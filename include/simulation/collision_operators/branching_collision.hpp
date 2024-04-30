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
#ifndef BRANCHING_COLLISION_H
#define BRANCHING_COLLISION_H

#include <materials/material.hpp>
#include <simulation/collision_operators/collision_operator.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/russian_roulette.hpp>

template <class FissOp>
class BranchingCollision {
 public:
  BranchingCollision() = default;

  void write_output_info(const std::string& base) const {
    if (mpi::rank != 0) return;

    auto& h5 = Output::instance().h5();
    h5.createAttribute<std::string>(base + "collision-operator", "branching");
  }

  void collision(Particle& p, MaterialHelper& mat,
                 ThreadLocalScores& thread_scores) const {
    // Sample a nuclide for the collision
    std::pair<const Nuclide*, MicroXSs> nuclide_info =
        mat.sample_nuclide(p.E(), p.rng);
    const auto& nuclide = nuclide_info.first;
    const auto& xs = nuclide_info.second;

    // Score kabs
    double k_abs_scr = p.wgt() * xs.nu_total * xs.fission / xs.total;
    thread_scores.k_abs_score += k_abs_scr;

    // Make all fission neutrons
    fission_operator.fission(p, xs, *nuclide);

    // Implicit capture
    p.set_weight(p.wgt() * (1. - (xs.absorption + xs.noise_copy) / xs.total));
    p.set_weight2(p.wgt2() * (1. - (xs.absorption + xs.noise_copy) / xs.total));

    // Roulette
    russian_roulette(p);

    // Scatter particle
    if (p.is_alive()) {
      // Perform scatter with nuclide
      scatter(p, *nuclide, xs);

      if (p.E() < settings::min_energy) {
        p.kill();
      }
    }
  }

 private:
  FissOp fission_operator;

  void scatter(Particle& p, const Nuclide& nuclide, const MicroXSs& xs) const {
    ScatterInfo sinfo = nuclide.sample_scatter(p.E(), p.u(), xs, p.rng);

    if (sinfo.yield == 0.) {
      // this scatter had a yield of 0... we have no choice but to kill
      // the particle now.
      p.kill();
      return;
    }

    // If yield is an integral quantity, we produce secondary neutrons in the
    // secondary bank, to keep weights from being too high and causing
    // variance issues.
    if (std::floor(sinfo.yield) == sinfo.yield && sinfo.yield != 1.) {
      int n = static_cast<int>(sinfo.yield) - 1;
      for (int i = 0; i < n; i++) {
        // Sample outgoing info
        ScatterInfo ninfo = nuclide.sample_scatter_mt(sinfo.mt, p.E(), p.u(),
                                                      xs.energy_index, p.rng);

        p.make_secondary(ninfo.direction, ninfo.energy,
                         p.wgt() * ninfo.weight_modifier,
                         p.wgt2() * ninfo.weight_modifier);
      }

      // This is set to 1 so that we have the correct weight for the last
      // neutron that we make outside this if block
      sinfo.yield = 1.;
    }

    p.set_direction(sinfo.direction);
    p.set_energy(sinfo.energy);
    p.set_weight(p.wgt() * sinfo.yield * sinfo.weight_modifier);
    p.set_weight2(p.wgt2() * sinfo.yield * sinfo.weight_modifier);
  }
};

#endif
