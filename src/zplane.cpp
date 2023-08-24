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
#include <geometry/surfaces/zplane.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

ZPlane::ZPlane(double z, BoundaryType bound, uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, z0{z} {}

int ZPlane::sign(const Position& r, const Direction& u) const {
  if (r.z() - z0 > SURFACE_COINCIDENT)
    return 1;
  else if (r.z() - z0 < -SURFACE_COINCIDENT)
    return -1;
  else {
    if (u.dot(norm(r)) > 0.)
      return 1;
    else
      return -1;
  }
}

double ZPlane::distance(const Position& r, const Direction& u,
                        bool on_surf) const {
  const double diff = z0 - r.z();
  if (on_surf || std::abs(diff) < SURFACE_COINCIDENT || u.z() == 0.)
    return INF;
  else if (diff / u.z() < 0.)
    return INF;
  else
    return diff / u.z();
}

Direction ZPlane::norm(const Position& /*r*/) const { return {0., 0., 1.}; }

//===========================================================================
// Non-Member Functions
std::shared_ptr<ZPlane> make_zplane(const YAML::Node& surface_node) {
  // Variables for surface
  double z0 = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get z0
  if (surface_node["z0"])
    z0 = surface_node["z0"].as<double>();
  else {
    fatal_error("ZPlane surface must have z0 defined.");
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
  return std::make_shared<ZPlane>(z0, boundary, id, name);
}
