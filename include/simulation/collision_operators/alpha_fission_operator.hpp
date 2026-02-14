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
#ifndef ALPHA_FISSION_OPERATOR_H
#define ALPHA_FISSION_OPERATOR_H

#include <materials/nuclide.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <cmath>

class AlphaFissionOperator {
 public:
  AlphaFissionOperator() = default;

  void fission(Particle& p, const MicroXSs& xs, const Nuclide& nuc,
               double alpha) const {
    const int n_new = this->n_fission_neutrons(p, xs);

    // Probability of a fission neutron being a delayed neutron.
    const double P_delayed = xs.nu_delayed / xs.nu_total;

    for (int i = 0; i < n_new; i++) {
      auto finfo =
          nuc.sample_fission(p.E(), p.u(), xs.energy_index, P_delayed, p.rng);

      // If parent was an alpha particle, we start with their weight
      double wgt = p.wgt();

      if (finfo.delayed) {
        const double lambda = finfo.precursor_decay_constant;
        wgt *= lambda / (lambda + alpha);
      }

      // Construct the new fission particle
      BankedParticle fiss_particle{p.r(),
                                   finfo.direction,
                                   finfo.energy,
                                   wgt,
                                   p.wgt2(),
                                   p.history_id(),
                                   p.daughter_counter(),
                                   p.family_id(),
                                   p.previous_collision_virtual(),
                                   p.previous_r(),
                                   p.previous_u(),
                                   p.previous_E(),
                                   p.E(),
                                   p.Esmp()};

      this->save_fission_particle(p, fiss_particle);
    }
  }

 private:
  int n_fission_neutrons(Particle& p, const MicroXSs& xs) const {
    // Since we are American, we cannot scale by keff here. Neutron production
    // is already scaled by keff
    return static_cast<int>(
        std::floor(xs.nu_total * xs.fission / xs.total + p.rng()));
  }

  void save_fission_particle(Particle& p, BankedParticle& fp) const {
    p.add_fission_particle(fp);
  }
};

#endif
