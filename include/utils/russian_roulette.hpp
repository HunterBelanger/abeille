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
#ifndef RUSSIAN_ROULETTE_H
#define RUSSIAN_ROULETTE_H

#include <simulation/particle.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <cmath>

void russian_roulette(Particle& p) {
  // Roulette first weight
  if (std::abs(p.wgt()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt()) / settings::wgt_survival);
    if (p.rng() < P_kill)
      p.set_weight(0.);
    else {
      p.set_weight(std::copysign(settings::wgt_survival, p.wgt()));
    }
  }

  // Roulette second weight
  if (std::abs(p.wgt2()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt2()) / settings::wgt_survival);
    if (p.rng() < P_kill)
      p.set_weight2(0.);
    else {
      p.set_weight2(std::copysign(settings::wgt_survival, p.wgt2()));
    }
  }

  // Check if we are still alive
  if (p.wgt() == 0. && p.wgt2() == 0.) p.kill();
}

#endif
