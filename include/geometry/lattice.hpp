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
#ifndef LATTICE_H
#define LATTICE_H

#include <cstdint>
#include <geometry/universe.hpp>

#include <yaml-cpp/yaml.h>

#include <map>
#include <vector>

class Lattice;

//===========================================================================
// Externals from geometry
namespace geometry {
extern std::vector<std::shared_ptr<Universe>> universes;
}  // namespace geometry

//===========================================================================
// Externals from parser
extern std::map<uint32_t, size_t> universe_id_to_indx;

extern void find_universe(const YAML::Node& input, uint32_t id);

class Lattice : public Universe {
 public:
  Lattice(uint32_t i_id, std::string i_name);
  virtual ~Lattice() = default;

  // Returns ture if position is inside a true lattice element
  // (not nullptr element), and false if it is outside or in a
  // dummy lattice element
  virtual bool is_inside(Position r, Direction u) const = 0;

  // Finds lattice element containing given position, transforms
  // coordinated to that element's frame, then asks that universe
  // for the cell of the local coordiante given.
  virtual UniqueCell get_cell(Position r, Direction u,
                              int32_t on_surf) const = 0;

  virtual UniqueCell get_cell(std::vector<GeoLilyPad>& stack, Position r,
                              Direction u, int32_t on_surf) const = 0;

  virtual void set_elements(std::vector<int32_t> univs) = 0;

  void set_outisde_universe(int32_t univ);

  Universe* outer_universe() const;

  std::size_t size() const { return lattice_universes.size(); }

  const Universe* get_universe(std::size_t ind) const;

  Universe* get_universe(std::size_t ind);

  Boundary get_boundary_condition(const Position& r, const Direction& u,
                                  int32_t on_surf) const override final;

  Boundary lost_get_boundary(const Position& r, const Direction& u,
                             int32_t on_surf) const override final;

  // get all material cells in universe
  std::set<uint32_t> get_all_mat_cells() const override final;

  // get number of cell instances across all universes
  uint32_t get_num_cell_instances(uint32_t cell_id) const override final;

  void make_offset_map() override final;

  bool contains_universe(uint32_t id) const override;

 protected:
  // Any lattice element that is negative gets directed to the
  // outer_universe. If outer_universe_index is also negative,
  // the particle is considered to be lost.
  std::vector<int32_t> lattice_universes;
  int32_t outer_universe_index;

};  // Lattice

//===========================================================================
// Non-Member functions
void make_lattice(const YAML::Node& uni_node, const YAML::Node& input);

#endif
