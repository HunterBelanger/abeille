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
#ifndef MODIFIED_FIXED_SOURCE_FISSION_OPERATOR_H
#define MODIFIED_FIXED_SOURCE_FISSION_OPERATOR_H

#include <materials/nuclide.hpp>
#include <simulation/collision_operators/fission_bank_saver.hpp>
#include <simulation/collision_operators/fission_operator.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/rng.hpp>

#include <cmath>
#include <memory>

class ModifiedFixedSourceFissionOperator
    : public FissionOperator<ModifiedFixedSourceFissionOperator,
                             FissionBankSaver> {
 public:
  ModifiedFixedSourceFissionOperator()
      : FissionOperator<ModifiedFixedSourceFissionOperator, FissionBankSaver>(
            FissionBankSaver()) {}

  int n_fission_neutrons(const Particle& p, const MicroXSs& xs,
                         pcg32& rng) const {
    const double k_abs_scr = p.wgt() * xs.nu_total * xs.fission / xs.total;

    // In fixed-source problems, we don't normalize the number of fission
    // particles produced. Problem must be sub critical !
    return static_cast<int>(std::floor(std::abs(k_abs_scr) + RNG::rand(rng)));
  }
};

#endif