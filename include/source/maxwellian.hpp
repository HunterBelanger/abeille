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
#ifndef MAXWELLIAN_H
#define MAXWELLIAN_H

#include <source/energy_distribution.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class Maxwellian : public EnergyDistribution {
 public:
  Maxwellian(double a);

  double sample(RNG& rng) const override final;

 private:
  double a_;  // Temperature in MeV
};

std::shared_ptr<Maxwellian> make_maxwellian(const YAML::Node& node);

#endif
