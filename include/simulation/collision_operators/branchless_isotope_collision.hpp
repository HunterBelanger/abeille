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
#ifndef BRANCHLESS_ISOTOPE_COLLISION_H
#define BRANCHLESS_ISOTOPE_COLLISION_H

#include <materials/material.hpp>
#include <simulation/collision_operators/collision_operator.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/russian_roulette.hpp>

template <class FissionSaver>
class BranchlessIsotopeCollision {
 public:
  BranchlessIsotopeCollision() = default;

  bool splitting() const { return splitting_; }
  void set_splitting(bool splt) { splitting_ = splt; }

  void write_output_info(const std::string& base) const {
    if (mpi::rank != 0) return;

    auto& h5 = Output::instance().h5();
    h5.createAttribute<std::string>(base + "collision-operator",
                                    "branchless-isotope");
    h5.createAttribute(base + "branchless-isotope-splitting", splitting_);
  }

  void collision(Particle& p, MaterialHelper& mat,
                 ThreadLocalScores& thread_scores) const {
    // Sample a nuclide for the collision
    std::pair<const Nuclide*, MicroXSs> nuclide_info =
        mat.sample_nuclide(p.E(), p.rng);
    const Nuclide& nuclide = *nuclide_info.first;
    const MicroXSs& microxs = nuclide_info.second;

    // Score contribution to kabs
    const double k_abs_scr =
        p.wgt() * microxs.nu_total * microxs.fission / microxs.total;
    thread_scores.k_abs_score += k_abs_scr;

    // Calcuate scatter probability and weight multiplier
    const double Pscatter = (microxs.total - microxs.absorption) /
                            (microxs.nu_total * microxs.fission +
                             (microxs.total - microxs.absorption));
    const double m = (microxs.nu_total * microxs.fission +
                      (microxs.total - microxs.absorption)) /
                     microxs.total;

    // Roulette
    russian_roulette(p);

    if (p.is_alive()) {
      if (p.rng() < Pscatter) {
        // Do scatter
        ScatterInfo sinfo =
            nuclide.sample_scatter(p.E(), p.u(), microxs, p.rng);

        if (sinfo.yield == 0.) {
          // this scatter had a yield of 0... we have no choice but to kill
          // the particle now.
          p.kill();
          return;
        }

        p.set_energy(sinfo.energy);
        p.set_direction(sinfo.direction);
        p.set_weight(p.wgt() * m * sinfo.yield);
        p.set_weight2(p.wgt2() * m * sinfo.yield);

        // Kill if energy is too low
        if (p.E() < settings::min_energy) {
          p.kill();
        }

        // Split particle if weight magnitude is too large
        if (splitting_ && p.is_alive() &&
            std::abs(p.wgt()) >= settings::wgt_split) {
          int n_new = static_cast<int>(std::floor(std::abs(p.wgt())));
          p.split(n_new);
        }
      } else {
        // Do fission
        double P_delayed = microxs.nu_delayed / microxs.nu_total;
        FissionInfo finfo = nuclide.sample_fission(
            p.E(), p.u(), microxs.energy_index, P_delayed, p.rng);

        // Construct the new fission particle
        BankedParticle fiss_particle{p.r(),
                                     finfo.direction,
                                     finfo.energy,
                                     p.wgt() * m,
                                     p.wgt2() * m,
                                     p.history_id(),
                                     p.daughter_counter(),
                                     p.family_id(),
                                     p.previous_collision_virtual(),
                                     p.previous_r(),
                                     p.previous_u(),
                                     p.previous_E(),
                                     p.E(),
                                     p.Esmp()};
        fiss_saver.save_fission_particle(p, fiss_particle);
        p.kill();
      }
    }
  }

 private:
  bool splitting_ = false;
  FissionSaver fiss_saver;
};

#endif
