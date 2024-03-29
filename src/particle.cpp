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
#include <simulation/particle.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <sstream>

Particle::Particle(Position r, Direction u, double engy, double wgt,
                   uint64_t id)
    : rng(),
      state{r, u, engy, wgt, 0.},
      history_id_(id),
      secondaries(),
      history_fission_bank(),
      r_birth_(state.position) {
  if (std::isnan(wgt)) {
    std::stringstream mssg;
    mssg << "Particle with id = " << id << " was instantiated with ";
    mssg << "wgt = NaN.";
    fatal_error(mssg.str());
  }
}

Particle::Particle(Position r, Direction u, double engy, double wgt,
                   double wgt2, uint64_t id)
    : rng(),
      state{r, u, engy, wgt, wgt2},
      history_id_(id),
      secondaries(),
      history_fission_bank(),
      r_birth_(state.position) {
  if (std::isnan(wgt)) {
    std::stringstream mssg;
    mssg << "Particle with id = " << id << " was instantiated with ";
    mssg << "wgt = NaN.";
    fatal_error(mssg.str());
  }

  if (std::isnan(wgt2)) {
    std::stringstream mssg;
    mssg << "Particle with id = " << id << " was instantiated with ";
    mssg << "wgt2 = NaN.";
    fatal_error(mssg.str());
  }
}
