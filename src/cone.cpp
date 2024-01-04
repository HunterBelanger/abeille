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
#include <source/cone.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <memory>

Cone::Cone(Direction u, double aperture)
    : direction_(u), aperture_(std::cos(aperture)) {}

Direction Cone::sample(pcg32& rng) const {
  // Sample mu in [aperture_, 1]
  double mu = (1. - aperture_) * RNG::rand(rng) + aperture_;

  // Sample phi in [0, 2pi]
  double phi = 2. * PI * RNG::rand(rng);

  return rotate_direction(direction_, mu, phi);
}

double Cone::aperture() const { return std::acos(aperture_); }

std::shared_ptr<Cone> make_cone_distribution(const YAML::Node& node) {
  if (!node["direction"] || !node["direction"].IsSequence() ||
      !(node["direction"].size() == 3)) {
    fatal_error("No valid direction entry for cone distribution.");
  }

  double x = node["direction"][0].as<double>();
  double y = node["direction"][1].as<double>();
  double z = node["direction"][2].as<double>();

  Direction u(x, y, z);

  if (!node["aperture"] || !node["aperture"].IsScalar()) {
    fatal_error("No valid aperture entry for cone distribution.");
  }
  double aperture = node["aperture"].as<double>();

  return std::make_shared<Cone>(u, aperture);
}
