/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#include <materials/material.hpp>
#include <materials/material_helper.hpp>
#include <simulation/tracker.hpp>
#include <simulation/transport_operators/carter_tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/majorant.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/energy_grid.hpp>

#include <xtensor/xtensor.hpp>

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

CarterTracker::CarterTracker(const std::vector<double>& mgxs)
    : EGrid(nullptr), Esmp(nullptr) {
  if (settings::energy_mode == settings::EnergyMode::CE) {
    fatal_error(
        "Using multi-group constructor for continuous energy simulation.");
  }

  if (mgxs.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Number of provided sampling cross section values does not agree "
            "with the number of energy groups.";
    fatal_error(mssg.str());
  }

  // Make the energy grid vector from the group bounds vector, and fill xs
  std::vector<double> egrid, xs;
  egrid.reserve(2 * settings::ngroups);
  xs.reserve(2 * settings::ngroups);
  std::size_t grp = 0;
  for (std::size_t i = 0; i < settings::energy_bounds.size(); i++) {
    egrid.push_back(settings::energy_bounds[i]);

    if (i != 0 && i != (settings::energy_bounds.size() - 1)) {
      egrid.push_back(settings::energy_bounds[i]);
    }

    if (grp < settings::ngroups) {
      const auto xsval = mgxs[grp];

      if (xsval <= 0.) {
        fatal_error("Sampling cross section must be > 0.");
      }

      xs.push_back(xsval);
      xs.push_back(xsval);
      grp++;
    }
  }

  EGrid = std::make_shared<pndl::EnergyGrid>(egrid);
  Esmp = std::make_shared<pndl::CrossSection>(xs, EGrid, 0);
}

CarterTracker::CarterTracker(const std::vector<double>& egrid,
                             const std::vector<double>& xs)
    : EGrid(nullptr), Esmp(nullptr) {
  if (settings::energy_mode == settings::EnergyMode::MG) {
    fatal_error(
        "Using continuous energy constructor for multi-group simulation.");
  }

  if (xs.size() != egrid.size()) {
    std::stringstream mssg;
    mssg << "Size of energy grid does not agree with size of cross section "
            "array.";
    fatal_error(mssg.str());
  }

  // Make sure energy grid is sorted
  if (std::is_sorted(egrid.begin(), egrid.end()) == false) {
    fatal_error("Sampling cross section energy grid is not sorted.");
  }

  // Make sure sampling xs is positive
  for (const auto& v : xs) {
    if (v <= 0.) {
      fatal_error("Sampling cross section must be > 0.");
    }
  }

  EGrid = std::make_shared<pndl::EnergyGrid>(egrid);
  Esmp = std::make_shared<pndl::CrossSection>(xs, EGrid, 0);
}

void CarterTracker::write_output_info(const std::string& base) const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();
  h5.createAttribute<std::string>(base + "transport-operator",
                                  "carter-tracking");

  // We now create a temporary array, which will hold the sampling xs info,
  // so we can save it to the output file.
  xt::xtensor<double, 2> smp_xs;
  smp_xs.resize({2, EGrid->size()});
  smp_xs.fill(0.);
  for (std::size_t i = 0; i < EGrid->size(); i++) {
    smp_xs(0, i) = (*EGrid)[i];
    smp_xs(1, i) = (*Esmp)[i];
  }
  std::vector<std::size_t> shape(smp_xs.shape().begin(), smp_xs.shape().end());
  auto smp_xs_ds =
      h5.createDataSet<double>(base + "sampling-xs", H5::DataSpace(shape));
  smp_xs_ds.write_raw(smp_xs.data());
}

void CarterTracker::transport(Particle& p, Tracker& trkr, MaterialHelper& mat,
                              ThreadLocalScores& thread_scores) const {
  bool had_collision = false;
  while (p.is_alive() && had_collision == false) {
    bool crossed_boundary = false;
    auto maj_indx = EGrid->get_lower_index(p.E());
    double Esample = Esmp->evaluate(p.E(), maj_indx) + mat.Ew(p.E());
    p.set_Esmp(Esample);  // Sampling XS saved for cancellation
    double d_coll = p.rng.exponential(Esample);
    Boundary bound(INF, -1, BoundaryType::Normal);

    // Score track length tally for boundary distance.
    // This is here because flux-like tallies are allowed with DT.
    // No other quantity should be scored with a TLE, as an error
    // should have been thrown when building all tallies.
    Tallies::instance().score_flight(p, trkr, std::min(d_coll, bound.distance),
                                     mat);

    // Try moving the distance to collision, and see if we land in a valid
    // material.
    trkr.move(d_coll);
    trkr.get_current();

    if (trkr.is_lost()) {
      // We got lost. This means we probably flew through a boundary
      // condition. We now go back, and actually find the B.C.
      trkr.set_r(p.r());
      trkr.restart_get_current();
      bound = trkr.get_boundary_condition();
      crossed_boundary = true;
    }

    if (crossed_boundary) {
      if (bound.boundary_type == BoundaryType::Vacuum) {
        p.kill();
        thread_scores.leakage_score += p.wgt();
        Position r_leak = p.r() + bound.distance * p.u();
        thread_scores.mig_score +=
            p.wgt() * (r_leak - p.r_birth()) * (r_leak - p.r_birth());
      } else if (bound.boundary_type == BoundaryType::Reflective) {
        trkr.do_reflection(p, bound);
        // Check if we are lost
        if (trkr.is_lost()) {
          std::stringstream mssg;
          mssg << "Particle " << p.history_id() << ".";
          mssg << p.secondary_id() << " has become lost.\n";
          mssg << "Previous valid coordinates: r = " << p.previous_r();
          mssg << ", u = " << p.previous_u() << ".\n";
          mssg << "Attempted reflection with surface "
               << geometry::surfaces[static_cast<std::size_t>(
                                         bound.surface_index)]
                      ->id();
          mssg << " at a distance of " << bound.distance << " cm.\n";
          mssg << "Currently lost at r = " << trkr.r() << ", u = " << trkr.u()
               << ".";
          fatal_error(mssg.str());
        }
        bound = trkr.get_boundary_condition();
      } else {
        fatal_error("Help me, how did I get here ?");
      }
    } else {
      // Update Position
      p.move(d_coll);
      mat.set_material(trkr.material(), p.E());

      // Get true cross section here
      const double Et = mat.Et(p.E());

      if (Esample >= Et) {
        const double Preal = Et / Esample;
        if (p.rng() < Preal) {
          // Flag real collision
          had_collision = true;
        }
      } else {
        const double D = Et / (2. * Et - Esample);
        const double F = Et / (D * Esample);

        if (p.rng() < D) {
          p.set_weight(p.wgt() * F);
          had_collision = true;
        } else {
          p.set_weight(-p.wgt() * F);
        }
      }
    }

    if (p.is_alive() && had_collision == false) {
      // Virtual collision
      p.set_previous_collision_virtual();
    }

    if (p.is_alive() && std::abs(p.wgt()) >= settings::wgt_split) {
      // Split particle, weight magnitude is too large
      int n_new = static_cast<int>(std::ceil(std::abs(p.wgt())));
      p.split(n_new);
    }
  }
}
