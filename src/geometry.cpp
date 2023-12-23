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
#include <geometry/geometry.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <memory>

namespace geometry {

//==========================================================================
// Global Variable Declrations
std::vector<std::shared_ptr<Surface>> surfaces;
std::vector<std::shared_ptr<Cell>> cells;
std::shared_ptr<Universe> root_universe;
std::vector<std::shared_ptr<Universe>> universes;

//==========================================================================
// Function Definitions
UniqueCell get_cell(const Position& r, const Direction& u, int32_t on_surf) {
  // Ask root_universe for cell. If no cell is found, answer
  // will be a nullptr
  return root_universe->get_cell(r, u, on_surf);
}

UniqueCell get_cell(std::vector<GeoLilyPad>& stack, const Position& r,
                    const Direction& u, int32_t on_surf) {
  // Ask root_universe for cell. If no cell is found, answer
  // will be a nullptr
  return root_universe->get_cell(stack, r, u, on_surf);
}

int32_t id_to_token(int32_t id) {
  if (surface_id_to_indx.find(static_cast<uint32_t>(std::abs(id))) ==
      surface_id_to_indx.end())
    return 0;

  int32_t token = static_cast<int32_t>(
      surface_id_to_indx[static_cast<uint32_t>(std::abs(id))]);
  token += 1;

  return token;
}

void do_reflection(Particle& p, Boundary boundary) {
  // First, we should get the surface from the index
  if (boundary.surface_index < 0) {
    fatal_error("Surface index is less than zero in reflection.");
  }
  const std::shared_ptr<Surface>& surface =
      surfaces[static_cast<std::size_t>(boundary.surface_index)];

  // Get new Position object to temporarily contain the current position
  Position r_pre_refs = p.r();

  // In the event a reflection occurs imediately after another, we must
  // ensure that we use r to be the previous position, or the position before
  // the first reflection, as that is our previous true position which
  // resulted from a sampled collision
  if (p.is_reflected()) r_pre_refs = p.previous_r();
  Direction u = p.u();

  // Get position of particle on surface
  Position r_on_surf = p.r() + boundary.distance * u;

  Direction n = surface->norm(r_on_surf);  // Get norm of surface

  // Get new direction after reflection
  Vector new_dir = u - 2. * (u * n) * n;  // Calc new direction
  Direction u_new = Direction{new_dir.x(), new_dir.y(), new_dir.z()};

  // Calc new previous position before reflections
  double d = boundary.distance + (p.r() - r_pre_refs).norm();
  Position r_prev = r_on_surf - d * u_new;

  // Update particle MUST BE DONE IN THIS ORDER
  // as set_position will also change previous_r, so we must set that second
  p.set_position(r_on_surf);
  p.set_previous_r(r_prev);
  p.set_direction(u_new);
  p.set_reflected(true);
}

}  // namespace geometry
