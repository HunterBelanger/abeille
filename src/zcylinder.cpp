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
#include <geometry/surfaces/zcylinder.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

ZCylinder::ZCylinder(double x_, double y_, double r_, BoundaryType bound,
                     uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, x0{x_}, y0{y_}, R{r_} {}

int ZCylinder::sign(const Position& r, const Direction& u) const {
  double x = r.x() - x0;
  double y = r.y() - y0;
  double eval = y * y + x * x - R * R;
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

double ZCylinder::distance(const Position& r, const Direction& u,
                           bool on_surf) const {
  double a = u.y() * u.y() + u.x() * u.x();
  if (a == 0.) return INF;

  double x = r.x() - x0;
  double y = r.y() - y0;
  double k = y * u.y() + x * u.x();
  double c = y * y + x * x - R * R;
  double quad = k * k - a * c;

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
    double d = (-k - std::sqrt(quad)) / a;
    if (d < 0.)
      return INF;
    else
      return d;
  }
}

Direction ZCylinder::norm(const Position& r) const {
  double x = r.x() - x0;
  double y = r.y() - y0;
  return {x, y, 0.};
}

//===========================================================================
// Non-Member Functions
std::shared_ptr<ZCylinder> make_zcylinder(const YAML::Node& surface_node) {
  // Variables for surface
  double x0 = 0., y0 = 0., r = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get x0
  if (surface_node["x0"])
    x0 = surface_node["x0"].as<double>();
  else {
    fatal_error("ZCylinder surface must have x0 defined.");
  }

  // Get y0
  if (surface_node["y0"])
    y0 = surface_node["y0"].as<double>();
  else {
    fatal_error("ZCylinder surface must have y0 defined.");
  }

  // Get r
  if (surface_node["r"])
    r = surface_node["r"].as<double>();
  else {
    fatal_error("ZCylinder surface must have r defined.");
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
  return std::make_shared<ZCylinder>(x0, y0, r, boundary, id, name);
}
