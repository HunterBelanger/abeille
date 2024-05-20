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
#include <geometry/surfaces/surface.hpp>
#include <materials/material_helper.hpp>
#include <simulation/tracker.hpp>
#include <simulation/transport_operators/surface_tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <sstream>

void SurfaceTracker::write_output_info(const std::string& base) const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();
  h5.createAttribute<std::string>(base + "transport-operator",
                                  "surface-tracking");
}

void SurfaceTracker::transport(Particle& p, Tracker& trkr, MaterialHelper& mat,
                               ThreadLocalScores& thread_scores) const {
  bool had_collision = false;

  while (p.is_alive() && had_collision == false) {
    double d_coll = p.rng.exponential(mat.Et(p.E()));
    auto bound = trkr.get_nearest_boundary();

    // score track length tally for boundary distance,
    // no matter what sort of boundary condition
    Tallies::instance().score_flight(p, trkr, std::min(d_coll, bound.distance),
                                     mat);

    double k_trk_scr =
        p.wgt() * std::min(d_coll, bound.distance) * mat.vEf(p.E());
    thread_scores.k_trk_score += k_trk_scr;

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
        trkr.cross_surface(bound);
        trkr.get_current();
        p.move(bound.distance);
        // Check if we are lost
        if (trkr.is_lost()) {
          std::stringstream mssg;
          mssg << "Particle " << p.history_id() << "." << p.secondary_id()
               << " has become lost.\n";
          mssg << " Previous valid coordinates: r = " << p.previous_r()
               << ", u = " << p.u() << ".\n";
          if (bound.surface_index >= 0) {
            mssg << " Attempted to cross surface "
                 << geometry::surfaces[static_cast<std::size_t>(
                                           bound.surface_index)]
                        ->id()
                 << " at a distance of " << bound.distance << " cm.\n";
          } else {
            mssg << " Attempted to cross surface at a distance of "
                 << bound.distance << " cm.\n";
          }
          mssg << " Currently lost at r = " << trkr.r() << ", u = " << trkr.u()
               << ".";
          fatal_error(mssg.str());
        }
        mat.set_material(trkr.material(), p.E());
      }
    } else {
      // Update Position
      p.move(d_coll);
      trkr.move(d_coll);

      // Flag real collision
      had_collision = true;
    }
  }
}
