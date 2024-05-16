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
#include <geometry/surfaces/all_surfaces.hpp>
#include <geometry/surfaces/surface.hpp>
#include <utils/error.hpp>

#include <string>

Surface::Surface(BoundaryType bound, uint32_t i_id, std::string i_name)
    : boundary_{bound}, id_{i_id}, name_{i_name} {}

BoundaryType Surface::boundary() const { return boundary_; }

uint32_t Surface::id() const { return id_; }

std::string Surface::name() const { return name_; }

std::shared_ptr<Surface> make_surface(const YAML::Node& surface_node) {
  // Try to get type
  std::string surf_type;
  if (surface_node["type"] && surface_node["type"].IsScalar())
    surf_type = surface_node["type"].as<std::string>();
  else {
    // Error, all surfaces must have a type
    fatal_error("Surface is missing \"type\" attribute.");
  }

  // Call appropriate function to build pointer to surface
  std::shared_ptr<Surface> surf_pntr = nullptr;
  if (surf_type == "xplane") {
    surf_pntr = make_xplane(surface_node);
  } else if (surf_type == "yplane") {
    surf_pntr = make_yplane(surface_node);
  } else if (surf_type == "zplane") {
    surf_pntr = make_zplane(surface_node);
  } else if (surf_type == "plane") {
    surf_pntr = make_plane(surface_node);
  } else if (surf_type == "xcylinder") {
    surf_pntr = make_xcylinder(surface_node);
  } else if (surf_type == "ycylinder") {
    surf_pntr = make_ycylinder(surface_node);
  } else if (surf_type == "zcylinder") {
    surf_pntr = make_zcylinder(surface_node);
  } else if (surf_type == "cylinder") {
    surf_pntr = make_cylinder(surface_node);
  } else if (surf_type == "sphere") {
    surf_pntr = make_sphere(surface_node);
  } else if (surf_type == "cross") {
    surf_pntr = make_cross(surface_node);
  } else {
    // Error, unknown surface type
    fatal_error("Surface type \"" + surf_type + "\" is unknown.");
  }

  return surf_pntr;
}