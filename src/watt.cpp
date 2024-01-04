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
#include <source/watt.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>

#include <cmath>

Watt::Watt(double a, double b) : EnergyDistribution(), max_(a), a_(a), b_(b) {
  if (a_ <= 0.) {
    fatal_error("Watt spectrum parameter a must be >= 0.");
  }

  if (b_ <= 0.) {
    fatal_error("Watt spectrum parameter b must be >= 0.");
  }
}

double Watt::sample(pcg32& rng) const {
  const double w = max_.sample(rng);

  return w + 0.25 * a_ * a_ * b_ +
         (2. * RNG::rand(rng) - 1.) * std::sqrt(a_ * a_ * b_ * w);
}

std::shared_ptr<Watt> make_watt(const YAML::Node& node) {
  if (!node["a"] || !node["a"].IsScalar()) {
    fatal_error("No valid \"a\" entry in watt distribution.");
  }
  const double a = node["a"].as<double>();
  if (a <= 0.) {
    fatal_error("Watt parameter a must be >= 0.");
  }

  if (!node["b"] || !node["b"].IsScalar()) {
    fatal_error("No valid \"b\" entry in watt distribution.");
  }
  const double b = node["b"].as<double>();
  if (b <= 0.) {
    fatal_error("Watt parameter b must be >= 0.");
  }

  return std::make_shared<Watt>(a, b);
}