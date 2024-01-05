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
#ifndef PARTICLE_MOVER_H
#define PARTICLE_MOVER_H

#include <materials/material_helper.hpp>
#include <noise_source/noise_maker.hpp>
#include <simulation/collision_operators/collision_operator.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <simulation/transport_operators/transport_operator.hpp>
#include <tallies/tallies.hpp>
#include <utils/error.hpp>
#include <utils/noise_parameters.hpp>
#include <utils/settings.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <optional>
#include <sstream>
#include <vector>

class IParticleMover {
  virtual ~IParticleMover() = default;

  virtual std::vector<BankedParticle> transport(std::vector<Particle>& bank, std::optional<NoiseParameters> noise = std::nullopt, std::vector<BankedParticle>* noise_bank = nullptr, const NoiseMaker* = nullptr) = 0;
};

template <TransportOperator TransportOp, CollisionOperator CollisionOp>
class ParticleMover : public IParticleMover {
 public:
  ParticleMover(TransportOp t, CollisionOp c) : transport_operator_(t), collision_operator_(c) {}

  std::vector<BankedParticle> transport(std::vector<Particle>& bank, std::optional<NoiseParameters> noise = std::nullopt, std::vector<BankedParticle>* noise_bank = nullptr, const NoiseMaker* noise_maker = nullptr) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel
#endif
    {
      // Thread local storage
      ThreadLocalScores thread_scores;

// Transport all particles in for thread
#ifdef ABEILLE_USE_OMP
#pragma omp for schedule(dynamic)
#endif
      for (std::size_t n = 0; n < bank.size(); n++) {
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
        // a nullptr. We also set the URR random variable here, if using
        // ptables.
        MaterialHelper mat(trkr.material(), p.E(), noise);
        if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);

        // While the partilce is alive, we continually do flights and collisions
        while (p.is_alive()) {
          transport_operator_.transport(p, trkr, mat, thread_scores);

          if (p.is_alive()) {
            // Sample noise source from the noise maker
            if (noise_maker) {
              noise_maker->sample_noise_source(p, mat);
            }

            // Score flux collision estimator with Sigma_t
            Tallies::instance().score_collision(p, mat);

            // Contribute to keff collision estimator and migration area scores
            thread_scores.k_col_score +=
                p.wgt() * mat.vEf(p.E()) / mat.Et(p.E());
            const double mig_dist = (p.r() - p.r_birth()).norm();
            thread_scores.mig_score +=
                p.wgt() * mat.Ea(p.E()) / mat.Et(p.E()) * mig_dist * mig_dist;

            // Perform the collision
            collision_operator_.collision(p, mat, thread_scores);

            // Update direction of the Tracker instance
            trkr.set_u(p.u());

            // Since the energy changed in the collision, we also need to update
            // the random values for the URR ptables, if being used.
            if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);

            // Set previous collision to real (for cancellation)
            p.set_previous_collision_real();
          }

          if (p.alive() == false) {
            // Particle was killed, either by the flight or collision.
            // Attempt a resurection
            resurection_ritual(p, trkr, mat);
          }
        }  // While alive
      }    // For all particles

      // Send all thread local scores to tallies instance
      Tallies::instance().score_k_col(thread_scores.k_col_score);
      Tallies::instance().score_k_abs(thread_scores.k_abs_score);
      Tallies::instance().score_k_trk(thread_scores.k_trk_score);
      Tallies::instance().score_k_tot(thread_scores.k_tot_score);
      Tallies::instance().score_leak(thread_scores.leakage_score);
      Tallies::instance().score_mig_area(thread_scores.mig_score);
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

 private:
  TransportOp transport_operator_;
  CollisionOp collision_operator_;

  void resurection_ritual(Particle& p, Tracker& trkr,
                          MaterialHelper& mat) const {
    p.resurect();

    if (p.is_alive()) {
      p.set_previous_collision_real();
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
};

std::shared_ptr<IParticleMover> make_particle_mover(const YAML::Node& sim, settings::SimMode mode);

#endif