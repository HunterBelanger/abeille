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
#include <source/direction_distribution.hpp>
#include <source/isotropic.hpp>
#include <source/mono_directional.hpp>
#include <utils/error.hpp>

#include <memory>

std::shared_ptr<DirectionDistribution> make_direction_distribution(
    const YAML::Node& node) {
  // Get the type of the distribution
  if (!node["type"] || !node["type"].IsScalar()) {
    fatal_error("No valid type provided to direction distribution entry.");
  }

  std::string type = node["type"].as<std::string>();

  std::shared_ptr<DirectionDistribution> dist = nullptr;

  if (type == "mono-directional") {
    dist = make_mono_directional_distribution(node);
  } else if (type == "isotropic") {
    dist = std::make_shared<Isotropic>();
  } else if (type == "cone") {
    dist = make_cone_distribution(node);
  } else {
    fatal_error("Invalid direction distribution type " + type + ".");
  }

  return dist;
}
