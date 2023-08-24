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
#include <geometry/surfaces/cylinder.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <cmath>

Cylinder::Cylinder(double x_, double y_, double z_, double u_, double v_,
                   double w_, double r_, BoundaryType bound, uint32_t i_id,
                   std::string i_name)
    : Surface{bound, i_id, i_name},
      x0{x_},
      y0{y_},
      z0{z_},
      u0{u_},
      v0{v_},
      w0{w_},
      R{r_},
      alpha{0.},
      beta{0.},
      gamma{0.} {
  // Normalize direction
  const double mag = std::sqrt(u0 * u0 + v0 * v0 + w0 * w0);

  if (mag == 0.) {
    fatal_error("Direction vector of Cylinder can not have a mangitude of 0.");
  }

  u0 /= mag;
  v0 /= mag;
  w0 /= mag;

  // Define alpha, beta, and gamma
  alpha = 1. - u0 * u0;
  beta = 1. - v0 * v0;
  gamma = 1. - w0 * w0;
}

int Cylinder::sign(const Position& r, const Direction& u) const {
  const double x = r.x() - x0;
  const double y = r.y() - y0;
  const double z = r.z() - z0;
  const double eval = alpha * x * x + beta * y * y + gamma * z * z - R * R;
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

double Cylinder::distance(const Position& r, const Direction& u,
                          bool on_surf) const {
  const double a =
      alpha * u.x() * u.x() + beta * u.y() * u.y() + gamma * u.z() * u.z();
  if (a == 0.) return INF;

  const double x = r.x() - x0;
  const double y = r.y() - y0;
  const double z = r.z() - z0;
  const double k = alpha * x * u.x() + beta * y * u.y() + gamma * z * u.z();
  const double c = alpha * x * x + beta * y * y + gamma * z * z - R * R;
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

Direction Cylinder::norm(const Position& r) const {
  return {alpha * (r.x() - x0), beta * (r.y() - y0), gamma * (r.z() - z0)};
}

//===========================================================================
// Non-Member Functions
std::shared_ptr<Cylinder> make_cylinder(const YAML::Node& surface_node) {
  // Variables for surface
  double x0 = 0., y0 = 0., z0 = 0., u0 = 0., v0 = 0., w0 = 0., r = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get x0
  if (surface_node["x0"])
    x0 = surface_node["x0"].as<double>();
  else {
    fatal_error("Cylinder surface must have x0 defined.");
  }

  // Get y0
  if (surface_node["y0"])
    y0 = surface_node["y0"].as<double>();
  else {
    fatal_error("Cylinder surface must have y0 defined.");
  }

  // Get z0
  if (surface_node["z0"])
    z0 = surface_node["z0"].as<double>();
  else {
    fatal_error("Cylinder surface must have z0 defined.");
  }

  // Get u0
  if (surface_node["u0"])
    u0 = surface_node["u0"].as<double>();
  else {
    fatal_error("Cylinder surface must have u0 defined.");
  }

  // Get v0
  if (surface_node["v0"])
    v0 = surface_node["v0"].as<double>();
  else {
    fatal_error("Cylinder surface must have v0 defined.");
  }

  // Get w0
  if (surface_node["w0"])
    w0 = surface_node["w0"].as<double>();
  else {
    fatal_error("Cylinder surface must have w0 defined.");
  }

  // Get r
  if (surface_node["r"])
    r = surface_node["r"].as<double>();
  else {
    fatal_error("Cylinder surface must have r defined.");
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
  return std::make_shared<Cylinder>(x0, y0, z0, u0, v0, w0, r, boundary, id,
                                    name);
}
