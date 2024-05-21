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
#include <geometry/surfaces/ycylinder.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

YCylinder::YCylinder(double x_, double z_, double r_, BoundaryType bound,
                     uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, x0{x_}, z0{z_}, R{r_} {}

int YCylinder::sign(const Position& r, const Direction& u) const {
  const double x = r.x() - x0;
  const double z = r.z() - z0;
  const double eval = x * x + z * z - R * R;
  if (eval > SURFACE_COINCIDENT)
    return 1;
  else if (eval < -SURFACE_COINCIDENT)
    return -1;
  else {
    if (u.dot(norm(r)) > 0.)
      return 1;
    else
      return -1;
  }
}

double YCylinder::distance(const Position& r, const Direction& u,
                           bool on_surf) const {
  const double a = u.x() * u.x() + u.z() * u.z();
  if (a == 0.) return INF;

  const double x = r.x() - x0;
  const double z = r.z() - z0;
  const double k = x * u.x() + z * u.z();
  const double c = x * x + z * z - R * R;
  const double quad = k * k - a * c;

  if (quad < 0.)
    return INF;
  else if (on_surf || std::abs(c) < SURFACE_COINCIDENT) {
    if (k >= 0.)
      return INF;
    else
      return (-k + std::sqrt(quad)) / a;
  } else if (c < 0.) {
    return (-k + std::sqrt(quad)) / a;
  } else {
    const double d = (-k - std::sqrt(quad)) / a;
    if (d < 0.)
      return INF;
    else
      return d;
  }
}

Direction YCylinder::norm(const Position& r) const {
  return {r.x() - x0, 0., r.z() - z0};
}

//===========================================================================
// Non-Member Functions
std::shared_ptr<YCylinder> make_ycylinder(const YAML::Node& surface_node) {
  // Variables for surface
  double x0 = 0., z0 = 0., r = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";
  
  // Get id
  if (surface_node["id"])
    id = surface_node["id"].as<uint32_t>();
  else {
    fatal_error(
        "Surface must have an id attribute with a unique positive integer.");
  }

  // Get x0
  if (surface_node["x0"] && surface_node["x0"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "ycylinder with id " << id << " has invalid x0 entry.";
    fatal_error(mssg.str());
  } else if (surface_node["x0"]) {
    x0 = surface_node["x0"].as<double>();
  }

  // Get z0
  if (surface_node["z0"] && surface_node["z0"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "ycylinder with id " << id << " has invalid z0 entry.";
    fatal_error(mssg.str());
  } else if (surface_node["z0"]) {
    z0 = surface_node["z0"].as<double>();
  }

  // Get r
  if (surface_node["r"])
    r = surface_node["r"].as<double>();
  else {
    fatal_error("xcylinder surface must have r defined.");
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

  // Get name
  if (surface_node["name"])
    name = surface_node["name"].as<std::string>();
  else
    name = "";

  // Make and return surface
  return std::make_shared<YCylinder>(x0, z0, r, boundary, id, name);
}
