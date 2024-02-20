/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
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
#ifndef TABULATED_ENERGY_H
#define TABULATED_ENERGY_H

#include <source/energy_distribution.hpp>

#include <PapillonNDL/pctable.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <vector>

// Assumes linear interpolation

class TabulatedEnergy : public EnergyDistribution {
 public:
  TabulatedEnergy(const pndl::PCTable& table);

  double sample(RNG& rng) const override final;

 private:
  pndl::PCTable table_;
};

std::shared_ptr<TabulatedEnergy> make_tabulated_energy(const YAML::Node& node);

#endif
