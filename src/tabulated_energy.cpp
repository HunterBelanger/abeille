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
#include <source/tabulated_energy.hpp>
#include <utils/error.hpp>
#include <utils/pctable.hpp>

#include <PapillonNDL/pndl_exception.hpp>

#include <algorithm>

TabulatedEnergy::TabulatedEnergy(const pndl::PCTable& table)
    : EnergyDistribution(), table_(table) {
  if (table_.min_value() <= 0.) {
    fatal_error("All values in a tabulated energy distribution must be > 0.");
  }
}

double TabulatedEnergy::sample(pcg32& rng) const {
  return table_.sample_value(RNG::rand(rng));
}

std::shared_ptr<TabulatedEnergy> make_tabulated_energy(const YAML::Node& node) {
  pndl::PCTable table = make_pctable(node);

  return std::make_shared<TabulatedEnergy>(table);
}