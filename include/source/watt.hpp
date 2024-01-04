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
#ifndef WATT_H
#define WATT_H

#include <source/energy_distribution.hpp>
#include <source/maxwellian.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class Watt : public EnergyDistribution {
 public:
  Watt(double a, double b);

  double sample(pcg32& rng) const override final;

 private:
  Maxwellian max_;
  double a_;  // Temperature in MeV
  double b_;  // Inverse MeV [1/MeV]
};

std::shared_ptr<Watt> make_watt(const YAML::Node& node);

#endif