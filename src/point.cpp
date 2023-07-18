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
#include <simulation/point.hpp>
#include <utils/error.hpp>

std::shared_ptr<Point> make_point_distribution(const YAML::Node& node) {
  if (!node["position"] || !node["position"].IsSequence() ||
      !(node["position"].size() == 3)) {
    fatal_error("No valid position entry for point spatial distribution.");
  }

  double x = node["position"][0].as<double>();
  double y = node["position"][1].as<double>();
  double z = node["position"][2].as<double>();

  Position r(x, y, z);

  return std::make_shared<Point>(r);
}
