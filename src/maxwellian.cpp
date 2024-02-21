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
#include <source/maxwellian.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>

Maxwellian::Maxwellian(double a) : EnergyDistribution(), a_(a) {
  if (a_ <= 0.) {
    fatal_error("Maxwellian parameter a must be >= 0.");
  }
}

double Maxwellian::sample(RNG& rng) const {
  const double xi1 = rng();
  const double xi2 = rng();
  const double xi3 = rng();

  const double c = std::cos(PI * xi3 / 2.);

  return -a_ * (std::log(xi1) + std::log(xi2) * c * c);
}

std::shared_ptr<Maxwellian> make_maxwellian(const YAML::Node& node) {
  if (!node["a"] || !node["a"].IsScalar()) {
    fatal_error("No valid \"a\" entry in maxwellian distribution.");
  }

  const double a = node["a"].as<double>();

  if (a <= 0.) {
    fatal_error("Maxwellian parameter a must be >= 0.");
  }

  return std::make_shared<Maxwellian>(a);
}
