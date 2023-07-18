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
#include <geometry/surfaces/yplane.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

YPlane::YPlane(double y, BoundaryType bound, uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, y0{y} {}

int YPlane::sign(const Position& r, const Direction& u) const {
  if (r.y() - y0 > SURFACE_COINCIDENT)
    return 1;
  else if (r.y() - y0 < -SURFACE_COINCIDENT)
    return -1;
  else {
    if (u.dot(norm(r)) > 0.)
      return 1;
    else
      return -1;
  }
}

double YPlane::distance(const Position& r, const Direction& u,
                        bool on_surf) const {
  double diff = y0 - r.y();
  if (on_surf || std::abs(diff) < SURFACE_COINCIDENT || u.y() == 0.)
    return INF;
  else if (diff / u.y() < 0.)
    return INF;
  else
    return diff / u.y();
}

Direction YPlane::norm(const Position& /*r*/) const { return {0., 1., 0.}; }

//===========================================================================
// Non-Member Functions
std::shared_ptr<YPlane> make_yplane(const YAML::Node& surface_node) {
  // Variables for surface
  double y0 = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get y0
  if (surface_node["y0"])
    y0 = surface_node["y0"].as<double>();
  else {
    fatal_error("YPlane surface must have y0 defined.");
  }

  // Get boundary type
  if (surface_node["boundary"]) {
    std::string boundary_string = surface_node["boundary"].as<std::string>();
    if (boundary_string == "vacuum")
      boundary = BoundaryType::Vacuum;
    else if (boundary_string == "reflective")
      boundary = BoundaryType::Reflective;
    else if (boundary_string == "normal")
      boundary = BoundaryType::Normal;
    else {
      fatal_error("Unknown boundary type \"" + boundary_string + "\".");
    }
  } else {
    boundary = BoundaryType::Normal;
  }

  // Get id
  if (surface_node["id"])
    id = surface_node["id"].as<uint32_t>();
  else {
    fatal_error(
        "Surface must have an id attribute with a unique positive integer.");
  }

  // Get name
  if (surface_node["name"])
    name = surface_node["name"].as<std::string>();
  else
    name = "";

  // Make and return surface
  return std::make_shared<YPlane>(y0, boundary, id, name);
}
