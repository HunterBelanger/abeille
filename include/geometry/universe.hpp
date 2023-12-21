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
#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <geometry/boundary.hpp>
#include <geometry/cell.hpp>
#include <geometry/geo_lily_pad.hpp>

#include <cstdint>

class Universe {
 public:
  Universe(uint32_t i_id, std::string i_name);
  virtual ~Universe() = default;

  virtual UniqueCell get_cell(Position r, Direction u,
                              int32_t on_surf) const = 0;

  virtual UniqueCell get_cell(std::vector<GeoLilyPad>& stack, Position r,
                              Direction u, int32_t on_surf) const = 0;

  virtual Boundary lost_get_boundary(const Position& r, const Direction& u,
                                     int32_t on_surf) const = 0;

  virtual Boundary get_boundary_condition(const Position& r, const Direction& u,
                                          int32_t on_surf) const = 0;

  //============================================================================
  // Lattice Methods
  virtual std::array<int32_t, 3> get_tile(Position r, Direction u) const;

  // Given the position in the frame of the lattice (NOT THE FRAME OF THE
  // TILE!), the distance to the edge of the provided tile is returned.
  virtual double distance_to_tile_boundary(Position r_local, Direction u,
                                           std::array<int32_t, 3> tile) const;
  //============================================================================

  virtual std::set<uint32_t> get_all_mat_cells() const = 0;

  virtual uint32_t get_num_cell_instances(uint32_t cell_id) const = 0;

  virtual void make_offset_map() = 0;

  virtual bool contains_universe(uint32_t id) const = 0;

  bool has_boundary_conditions() const { return has_boundary_conditions_; }

  uint32_t id() const;

  std::string name() const;

 protected:
  std::vector<std::map<const uint32_t, uint32_t>> cell_offset_map;
  uint32_t id_;
  std::string name_;
  bool has_boundary_conditions_ = true;

};  // Universe

#endif
