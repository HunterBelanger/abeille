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
#include <simulation/transport_operators/implicit_leakage_delta_tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/majorant.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/energy_grid.hpp>

#include <ndarray.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

ImplicitLeakageDeltaTracker::ImplicitLeakageDeltaTracker()
    : EGrid(nullptr), Emaj(nullptr) {
  Output::instance().write(" Finding majorant cross sections.\n");
  auto Egrid_Emaj_pair = make_majorant_xs();
  EGrid = std::make_shared<pndl::EnergyGrid>(Egrid_Emaj_pair.first);
  Emaj = std::make_shared<pndl::CrossSection>(Egrid_Emaj_pair.second, EGrid, 0);
}

void ImplicitLeakageDeltaTracker::write_output_info(
    const std::string& base) const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();
  h5.createAttribute<std::string>(base + "transport-operator",
                                  "implicit-leakage-delta-tracking");

  // We now create a temporary array, which will hold the majorant xs info,
  // so we can save it to the output file.
  NDArray<double> maj_xs({2, EGrid->size()});
  for (std::size_t i = 0; i < EGrid->size(); i++) {
    maj_xs(0, i) = (*EGrid)[i];
    maj_xs(1, i) = (*Emaj)[i];
  }
  auto maj_xs_ds = h5.createDataSet<double>(base + "majorant-xs",
                                            H5::DataSpace(maj_xs.shape()));
  maj_xs_ds.write_raw(&maj_xs[0]);
}

void ImplicitLeakageDeltaTracker::transport(
    Particle& p, Tracker& trkr, MaterialHelper& mat,
    ThreadLocalScores& thread_scores) const {
  bool had_collision = false;
  while (p.is_alive() && had_collision == false) {
    auto maj_indx = EGrid->get_lower_index(p.E());
    double Emajorant = Emaj->evaluate(p.E(), maj_indx) + mat.Ew(p.E());
    p.set_Esmp(Emajorant);  // Sampling XS saved for cancellation
    auto bound = trkr.get_boundary_condition();

    double d_coll = 0.;
    if (bound.boundary_type == BoundaryType::Vacuum) {
      // Calculate probability of leaking and not leaking
      const double P_leak = std::exp(-Emajorant * bound.distance);
      const double P_no_leak = 1. - P_leak;

      // Get the weight which flys all the way to boundary and leaks
      const double wgt_leak = p.wgt() * P_leak;
      const double wgt2_leak = p.wgt2() * P_leak;

      // Get the weight which will have a collision
      const double wgt_collides = p.wgt() * P_no_leak;
      const double wgt2_collides = p.wgt2() * P_no_leak;

      // Score implicit leackage
      thread_scores.leakage_score += wgt_leak;

      // Score implicit migration area
      Position r_leak = p.r() + bound.distance * p.u();
      thread_scores.mig_score +=
          wgt_leak * (r_leak - p.r_birth()) * (r_leak - p.r_birth());

      // Sample flight distance
      d_coll = -std::log(1. - P_no_leak * p.rng()) / Emajorant;

      // Score the TLE for the portion which would have leaked
      p.set_weight(wgt_leak);
      p.set_weight2(wgt2_leak);
      Tallies::instance().score_flight(p, bound.distance, mat);
      Tallies::instance().score_flight(p, trkr, bound.distance,
                                     mat);

      // Reduce weight
      p.set_weight(wgt_collides);
      p.set_weight2(wgt2_collides);

      // Score TLE for the portion which only goes to the collision site
      Tallies::instance().score_flight(p, d_coll, mat);
      Tallies::instance().score_flight(p, trkr, d_coll,
                                     mat);
    } else {
      d_coll = p.rng.exponential(Emajorant);

      // Score track length tally for boundary distance.
      // This is here because flux-like tallies are allowed with DT.
      // No other quantity should be scored with a TLE, as an error
      // should have been thrown when building all tallies.
      Tallies::instance().score_flight(p, std::min(d_coll, bound.distance),
                                       mat);
      Tallies::instance().score_flight(p, trkr, std::min(d_coll, bound.distance),
                                     mat);
    }

    if (bound.distance < d_coll ||
        std::abs(bound.distance - d_coll) < BOUNDRY_TOL) {
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
      // Update Position
      p.move(d_coll);
      trkr.move(d_coll);
      trkr.get_current();
      // Check if we are lost
      if (trkr.is_lost()) {
        std::stringstream mssg;
        mssg << "Particle " << p.history_id() << ".";
        mssg << p.secondary_id() << " has become lost.\n";
        mssg << "Previous valid coordinates: r = " << p.previous_r();
        mssg << ", u = " << p.previous_u() << ".\n";
        mssg << "Attempted to fly a distance of " << d_coll << " cm.\n";
        mssg << "Currently lost at r = " << trkr.r() << ", u = " << trkr.u()
             << ".";
        fatal_error(mssg.str());
      }
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
