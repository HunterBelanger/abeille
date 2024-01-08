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
#ifndef KEFF_FISSION_OPERATOR_H
#define KEFF_FISSION_OPERATOR_H

#include <materials/nuclide.hpp>
#include <simulation/collision_operators/fission_bank_saver.hpp>
#include <simulation/collision_operators/fission_operator.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/rng.hpp>

#include <cmath>
#include <memory>

class KeffFissionOperator
    : public FissionOperator<KeffFissionOperator, FissionBankSaver> {
 public:
  KeffFissionOperator() = default;

  int n_fission_neutrons(Particle& p, const MicroXSs& xs) const {
    const double k_abs_scr = p.wgt() * xs.nu_total * xs.fission / xs.total;

    // In k-eigenvalue simulations, we normalize particle production by the
    // keff of the previous generation, so that the number of particles per
    // generation stays approximately constant.
    return static_cast<int>(std::floor(
        std::abs(k_abs_scr) / Tallies::instance().kcol() + RNG::rand(p.rng)));
  }
};

#endif