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
#ifndef CYLINDER_H
#define CYLINDER_H

#include <geometry/surfaces/surface.hpp>

#include <yaml-cpp/yaml.h>

class Cylinder : public Surface {
 public:
  Cylinder(double x_, double y_, double z_, double u_, double v_, double w_,
           double r_, BoundaryType bound, uint32_t i_id, std::string i_name);
  ~Cylinder() = default;

  int sign(const Position& r, const Direction& u) const override;

  double distance(const Position& r, const Direction& u,
                  bool on_surf) const override;

  Direction norm(const Position& r) const override;

 private:
  double x0, y0, z0, u0, v0, w0, R;
  double alpha, beta, gamma;

};  // Cylinder

//===========================================================================
// Non-Member Functions
std::shared_ptr<Cylinder> make_cylinder(const YAML::Node& surface_node);

#endif
