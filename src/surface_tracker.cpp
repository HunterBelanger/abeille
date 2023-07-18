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
#include <simulation/surface_tracker.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

#include <iostream>
#include <string>

#ifdef ABEILLE_USE_OMP
#include <omp.h>
#endif

std::vector<BankedParticle> SurfaceTracker::transport(
    std::vector<Particle>& bank, bool noise,
    std::vector<BankedParticle>* noise_bank, const NoiseMaker* noise_maker) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel
#endif
  {
    // Thread local storage
    ThreadLocalScores thread_scores;

// Transport all particles in for thread
#ifdef ABEILLE_USE_OMP
#pragma omp for schedule(static)
#endif
    for (size_t n = 0; n < bank.size(); n++) {
      // Particle and its personal tracker
      Particle& p = bank[n];
      Tracker trkr(p.r(), p.u());

      // If we got lost, kill the particle
      if (trkr.is_lost()) {
        std::stringstream mssg;
        mssg << "Particle become lost at " << p.r() << ", ";
        mssg << " u = " << p.u() << ", token = " << trkr.surface_token();
        warning(mssg.str());
        p.kill();
      }
      // Only make helper if we aren't lost, to make sure that material isn't
      // a nullptr. We also set the URR random variable here, if using ptables.
      MaterialHelper mat(trkr.material(), p.E());
      if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);

      while (p.is_alive()) {
        bool had_collision = false;
        double d_coll = RNG::exponential(p.rng, mat.Et(p.E(), noise));
        auto bound = trkr.get_nearest_boundary();

        // score track length tally for boundary distance,
        // no matter what sort of boundary condition
        tallies->score_flight(p, std::min(d_coll, bound.distance), mat,
                              settings::converged);
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
              mssg << "Currently lost at r = " << trkr.r()
                   << ", u = " << trkr.u() << ".";
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
              mssg << " Currently lost at r = " << trkr.r()
                   << ", u = " << trkr.u() << ".";
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

        if (p.is_alive() && had_collision) {  // real collision
          collision(p, mat, thread_scores, noise, noise_maker);
          trkr.set_u(p.u());
          if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);
        }  // If alive for real collision

        if (!p.is_alive()) {
          // Attempt a resurection
          p.resurect();

          if (p.is_alive()) {
            trkr.set_r(p.r());
            trkr.set_u(p.u());
            trkr.restart_get_current();
            // Check if we are lost
            if (trkr.is_lost()) {
              std::stringstream mssg;
              mssg << "Particle " << p.history_id() << ".";
              mssg << p.secondary_id() << " has become lost.\n";
              mssg << "Attempted resurection at r = " << trkr.r();
              mssg << ", u = " << trkr.u() << ".";
              fatal_error(mssg.str());
            }
            mat.set_material(trkr.material(), p.E());
            if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);
          } else if (settings::rng_stride_warnings) {
            // History is truly dead.
            // Check if we went past the particle stride.
            uint64_t n_rng_calls = p.number_of_rng_calls();
            if (n_rng_calls > settings::rng_stride) {
              // This isn't really a good thing. We should
              // write a warning.
              std::string mssg = "History " + std::to_string(p.history_id()) +
                                 " overran the RNG stride.";
              warning(mssg);
            }
          }
        }
      }  // While alive
    }    // For all particles

    // Send all thread local scores to tallies instance
    tallies->score_k_col(thread_scores.k_col_score);
    tallies->score_k_abs(thread_scores.k_abs_score);
    tallies->score_k_trk(thread_scores.k_trk_score);
    tallies->score_k_tot(thread_scores.k_tot_score);
    tallies->score_leak(thread_scores.leakage_score);
    tallies->score_mig_area(thread_scores.mig_score);
    thread_scores.k_col_score = 0.;
    thread_scores.k_abs_score = 0.;
    thread_scores.k_trk_score = 0.;
    thread_scores.k_tot_score = 0.;
    thread_scores.leakage_score = 0.;
    thread_scores.mig_score = 0.;
  }  // Parallel

  // Vector to contain all fission daughters for all threads
  std::vector<BankedParticle> fission_neutrons;

  // Empty all particle fission banks into the main one
  for (auto& p : bank) {
    p.empty_fission_bank(fission_neutrons);
  }

  if (noise_bank && noise_maker) {
    for (auto& p : bank) {
      p.empty_noise_bank(*noise_bank);
    }
  }

  // Can now clear the old bank
  bank.clear();

  return fission_neutrons;
}
