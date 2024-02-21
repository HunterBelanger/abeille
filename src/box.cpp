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
#include <source/box.hpp>
#include <utils/error.hpp>

#include <memory>

Box::Box(Position low, Position hi) : low_(low), hi_(hi) {
  if (low.x() > hi.x() || low.y() > hi.y() || low.z() > hi.z()) {
    fatal_error("Coordinate of low is greater than hi.");
  }
}

Position Box::sample(RNG& rng) const {
  double x = (hi_.x() - low_.x()) * rng() + low_.x();
  double y = (hi_.y() - low_.y()) * rng() + low_.y();
  double z = (hi_.z() - low_.z()) * rng() + low_.z();
  return {x, y, z};
}

std::shared_ptr<Box> make_box_distribution(const YAML::Node& node) {
  if (!node["low"] || !node["low"].IsSequence() || !(node["low"].size() == 3)) {
    fatal_error("No valid low entry for box spatial distribution.");
  }

  double xl = node["low"][0].as<double>();
  double yl = node["low"][1].as<double>();
  double zl = node["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  if (!node["hi"] || !node["hi"].IsSequence() || !(node["hi"].size() == 3)) {
    fatal_error("No valid hi entry for box spatial distribution.");
  }

  double xh = node["hi"][0].as<double>();
  double yh = node["hi"][1].as<double>();
  double zh = node["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  return std::make_shared<Box>(r_low, r_hi);
}
