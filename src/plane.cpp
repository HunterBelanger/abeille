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
#include <geometry/surfaces/plane.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

Plane::Plane(double A_, double B_, double C_, double D_, BoundaryType bound,
             uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, A{A_}, B{B_}, C{C_}, D{D_} {}

int Plane::sign(const Position& r, const Direction& u) const {
  double eval = A * r.x() + B * r.y() + C * r.z() - D;
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

double Plane::distance(const Position& r, const Direction& u,
                       bool on_surf) const {
  double num = D - A * r.x() - B * r.y() - C * r.z();
  double denom = A * u.x() + B * u.y() + C * u.z();
  double d = num / denom;
  if (on_surf || std::abs(d) < SURFACE_COINCIDENT || denom == 0.)
    return INF;
  else if (d < 0.)
    return INF;
  else
    return d;
}

Direction Plane::norm(const Position& /*r*/) const { return {A, B, C}; }

//===========================================================================
// Non-Member Functions
std::shared_ptr<Plane> make_plane(const YAML::Node& surface_node) {
  // Variables for surface
  double A = 0., B = 0., C = 0., D = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get A
  if (surface_node["A"])
    A = surface_node["A"].as<double>();
  else {
    fatal_error("Plane surface must have A defined.");
  }

  // Get B
  if (surface_node["B"])
    B = surface_node["B"].as<double>();
  else {
    fatal_error("Plane surface must have B defined.");
  }

  // Get C
  if (surface_node["C"])
    C = surface_node["C"].as<double>();
  else {
    fatal_error("Plane surface must have C defined.");
  }

  // Get D
  if (surface_node["D"])
    D = surface_node["D"].as<double>();
  else {
    fatal_error("Plane surface must have D defined.");
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
  return std::make_shared<Plane>(A, B, C, D, boundary, id, name);
}
