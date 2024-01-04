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
#ifndef FISSION_OPERATOR_H
#define FISSION_OPERATOR_H

#include <materials/nuclide.hpp>
#include <simulation/particle.hpp>
#include <utils/rng.hpp>

template <class T, class FissionSaver>
class FissionOperator {
 public:
  FissionOperator(FissionSaver fiss_saver): fission_saver(fiss_saver) {}

  void fission(Particle& p, const MicroXSs& xs, const Nuclide& nuc) const {
    const int n_new =
        static_cast<const T*>(this)->n_fission_neutrons(p, xs, p.rng);

    // Probability of a fission neutron being a delayed neutron.
    const double P_delayed = xs.nu_delayed / xs.nu_total;

    for (int i = 0; i < n_new; i++) {
      auto finfo =
          nuc.sample_fission(p.E(), p.u(), xs.energy_index, P_delayed, p.rng);

      // If we are doing transport of non-noise particles, the weight of fission
      // particles should be +/- 1, depending on sign of parent.
      double wgt = p.wgt() > 0. ? 1. : -1.;
      double wgt2 = 0.;

      // Construct the new fission particle
      BankedParticle fiss_particle{p.r(),
                                   finfo.direction,
                                   finfo.energy,
                                   wgt,
                                   wgt2,
                                   p.history_id(),
                                   p.daughter_counter(),
                                   p.family_id(),
                                   p.previous_collision_virtual(),
                                   p.previous_r(),
                                   p.previous_u(),
                                   p.previous_E(),
                                   p.E(),
                                   p.Esmp()};

      fission_saver.save_fission_particle(p, fiss_particle);
    }
  }

 private:
  FissionSaver fission_saver;
};

#endif