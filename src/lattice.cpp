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
#include <geometry/boundary.hpp>
#include <geometry/geometry.hpp>
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <utils/error.hpp>

#include <cstdint>

Lattice::Lattice(uint32_t i_id, std::string i_name)
    : Universe{i_id, i_name}, lattice_universes{}, outer_universe_index{-1} {
  this->has_boundary_conditions_ = false;

  if (this->outer_universe()) {
    this->has_boundary_conditions_ =
        this->outer_universe()->has_boundary_conditions();
  }
}

void Lattice::set_outisde_universe(int32_t univ) {
  outer_universe_index = univ;

  if (this->outer_universe()) {
    this->has_boundary_conditions_ =
        this->outer_universe()->has_boundary_conditions();
  }
}

Universe* Lattice::outer_universe() const {
  if (outer_universe_index < 0) return nullptr;

  return geometry::universes[static_cast<std::size_t>(outer_universe_index)]
      .get();
}

const Universe* Lattice::get_universe(std::size_t ind) const {
  if (lattice_universes[ind] < 0) return nullptr;

  std::size_t uni_indx = static_cast<std::size_t>(lattice_universes[ind]);

  return geometry::universes[uni_indx].get();
}

Universe* Lattice::get_universe(std::size_t ind) {
  if (lattice_universes[ind] < 0) return nullptr;

  std::size_t uni_indx = static_cast<std::size_t>(lattice_universes[ind]);

  return geometry::universes[uni_indx].get();
}

Boundary Lattice::get_boundary_condition(const Position& r, const Direction& u,
                                         int32_t on_surf) const {
  if (this->has_boundary_conditions()) {
    // Get lattice

    // If there is a boundary condition here, we know there is an outer
    // universe, so we are safe to dereference.
    return this->outer_universe()->get_boundary_condition(r, u, on_surf);
  }

  // No outer universe or no boundary conditions, so this goes out to infinity
  Boundary ret_bound(INF, -1, BoundaryType::Vacuum);
  ret_bound.token = 0;
  return ret_bound;
}

Boundary Lattice::lost_get_boundary(const Position& r, const Direction& u,
                                    int32_t on_surf) const {
  // Get lattice

  const bool is_inside = this->is_inside(r, u);
  if (this->outer_universe() && is_inside == false) {
    return this->outer_universe()->lost_get_boundary(r, u, on_surf);
  }

  std::array<int32_t, 3> tile = this->get_tile(r, u);
  Boundary ret_bound(this->distance_to_tile_boundary(r, u, tile), -1,
                     BoundaryType::Normal);
  ret_bound.token = 0;

  return ret_bound;
}

bool Lattice::contains_universe(uint32_t id) const {
  // Get lattice

  // First check outer
  if (this->outer_universe()) {
    if (this->outer_universe()->id() == id ||
        this->outer_universe()->contains_universe(id))
      return true;
  }

  // Now go through all other universes in lattice
  for (std::size_t ind = 0; ind < this->size(); ind++) {
    const Universe* uni = this->get_universe(ind);

    if (uni) {
      if (uni->id() == id || uni->contains_universe(id)) return true;
    }
  }

  return false;
}

uint32_t Lattice::get_num_cell_instances(uint32_t cell_id) const {
  uint32_t instances = 0;

  // go through all tiles/universes
  for (std::size_t ind = 0; ind < this->size(); ind++) {
    const Universe* uni = this->get_universe(ind);
    if (uni) {
      instances += uni->get_num_cell_instances(cell_id);
    }
  }

  if (this->outer_universe()) {
    instances += this->outer_universe()->get_num_cell_instances(cell_id);
  }

  return instances;
}

std::set<uint32_t> Lattice::get_all_mat_cells() const {
  std::set<uint32_t> mat_cells;

  // go through all tiles/universes
  for (std::size_t ind = 0; ind < this->size(); ind++) {
    const Universe* uni = this->get_universe(ind);
    if (uni) {
      auto uni_mat_cells = uni->get_all_mat_cells();
      mat_cells.insert(uni_mat_cells.begin(), uni_mat_cells.end());
    }
  }

  if (this->outer_universe()) {
    auto outer_mat_cells = this->outer_universe()->get_all_mat_cells();
    mat_cells.insert(outer_mat_cells.begin(), outer_mat_cells.end());
  }

  return mat_cells;
}

void Lattice::make_offset_map() {
  cell_offset_map.resize(this->size() + (this->outer_universe() ? 1 : 0));
  // Set of all material cells contained in this universe
  std::set<uint32_t> mat_cell_ids = this->get_all_mat_cells();

  // Set all offsets to be zero
  for (auto& map : cell_offset_map) {
    for (uint32_t mat_cell_id : mat_cell_ids) {
      map[mat_cell_id] = 0;
    }
  }

  // Now go through and build all offsets
  for (std::size_t i = 1; i < this->size(); i++) {
    const Universe* uni = this->get_universe(i - 1);

    for (uint32_t mat_cell_id : mat_cell_ids) {
      cell_offset_map[i][mat_cell_id] = cell_offset_map[i - 1][mat_cell_id];

      if (uni) {
        cell_offset_map[i][mat_cell_id] +=
            uni->get_num_cell_instances(mat_cell_id);
      }
    }
  }

  // Also get offsets for outer universe
  if (this->outer_universe()) {
    const std::size_t i = cell_offset_map.size() - 1;
    auto uni = this->get_universe(this->size() - 1);

    for (uint32_t mat_cell_id : mat_cell_ids) {
      cell_offset_map[i][mat_cell_id] = cell_offset_map[i - 1][mat_cell_id];

      if (uni) {
        cell_offset_map[i][mat_cell_id] +=
            uni->get_num_cell_instances(mat_cell_id);
      }
    }
  }
}

void make_lattice(const YAML::Node& uni_node, const YAML::Node& input) {
  if (!uni_node["type"] || !uni_node["type"].IsScalar()) {
    fatal_error("Lattice instances must have a valid type.");
  }

  std::string type = uni_node["type"].as<std::string>();

  if (type == "rectlinear") {
    make_rect_lattice(uni_node, input);
  } else if (type == "hexagonal") {
    make_hex_lattice(uni_node, input);
  } else {
    fatal_error("Unknown lattice type " + type + ".");
  }
}
