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
#include <simulation/transport_operators/delta_tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/majorant.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/energy_grid.hpp>

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

DeltaTracker::DeltaTracker() : EGrid(nullptr), Emaj(nullptr) {
  Output::instance().write(" Finding majorant cross sections.\n");
  auto Egrid_Emaj_pair = make_majorant_xs();
  EGrid = std::make_shared<pndl::EnergyGrid>(Egrid_Emaj_pair.first);
  Emaj = std::make_shared<pndl::CrossSection>(Egrid_Emaj_pair.second, EGrid, 0);
}

void DeltaTracker::write_output_info(const std::string& base) const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();
  h5.createAttribute<std::string>(base + "transport-operator",
                                  "delta-tracking");

  // We now create a temporary array, which will hold the majorant xs info,
  // so we can save it to the output file.
  xt::xtensor<double, 2> maj_xs;
  maj_xs.resize({2, EGrid->size()});
  maj_xs.fill(0.);
  for (std::size_t i = 0; i < EGrid->size(); i++) {
    maj_xs(0, i) = (*EGrid)[i];
    maj_xs(1, i) = (*Emaj)[i];
  }
  std::vector<std::size_t> shape(maj_xs.shape().begin(), maj_xs.shape().end());
  auto maj_xs_ds =
      h5.createDataSet<double>(base + "majorant-xs", H5::DataSpace(shape));
  maj_xs_ds.write_raw(maj_xs.data());
}

void DeltaTracker::transport(Particle& p, Tracker& trkr, MaterialHelper& mat,
                             ThreadLocalScores& thread_scores) const {
  bool had_collision = false;
  while (p.is_alive() && had_collision == false) {
    bool crossed_boundary = false;
    auto maj_indx = EGrid->get_lower_index(p.E());
    double Emajorant = Emaj->evaluate(p.E(), maj_indx) + mat.Ew(p.E());
    p.set_Esmp(Emajorant);  // Sampling XS saved for cancellation
    double d_coll = p.rng.exponential(Emajorant);
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
      } else {
        fatal_error("Help me, how did I get here ?");
      }
    } else {
      // Update position on the particle. Tracker is already up to date.
      p.move(d_coll);

      // Set the material for the current position.
      mat.set_material(trkr.material(), p.E());

      // Get true cross section here
      double Et = mat.Et(p.E());

      if (Et - Emajorant > 1.E-10) {
        std::stringstream mssg;
        mssg << "Total cross section excedeed majorant at ";
        mssg << p.E() << " MeV.";
        mssg << " Et = " << Et << ", Emaj = " << Emajorant << "\n";
        fatal_error(mssg.str());
      }

      if (p.rng() < (Et / Emajorant)) {
        // Flag real collision
        had_collision = true;
      }
    }

    if (p.is_alive() && had_collision == false) {
      // virtual collision
      p.set_previous_collision_virtual();
    }
  }
}
