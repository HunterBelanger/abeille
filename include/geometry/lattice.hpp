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
extern std::vector<std::shared_ptr<Lattice>> lattices;
extern std::vector<std::shared_ptr<Universe>> universes;
}  // namespace geometry

//===========================================================================
// Externals from parser
extern std::map<uint32_t, size_t> universe_id_to_indx;
extern std::map<uint32_t, size_t> lattice_id_to_indx;
extern void find_universe(const YAML::Node& input, uint32_t id);

class Lattice {
 public:
  Lattice(uint32_t i_id, std::string i_name);
  virtual ~Lattice() = default;

  // Returns ture if position is inside a true lattice element
  // (not nullptr element), and false if it is outside or in a
  // dummy lattice element
  virtual bool is_inside(Position r, Direction u) const = 0;

  virtual std::array<int32_t, 3> get_tile(Position r, Direction u) const = 0;

  // Finds lattice element containing given position, transforms
  // coordinated to that elements frame, then asks that universe
  // for the cell of the local coordiante given.
  virtual Cell* get_cell(Position r, Direction u, int32_t on_surf) const = 0;

  virtual Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r,
                         Direction u, int32_t on_surf) const = 0;

  // Given the position in the frame of the lattice (NOT THE FRAME OF THE
  // TILE!), the distance to the edge of the provided tile is returned.
  virtual double distance_to_tile_boundary(
      Position r_local, Direction u, std::array<int32_t, 3> tile) const = 0;

  virtual void set_elements(std::vector<int32_t> univs) = 0;

  void set_outisde_universe(int32_t univ);

  Universe* outer_universe() const;

  std::size_t size() const { return lattice_universes.size(); }

  const Universe* get_universe(std::size_t ind) const;

  uint32_t id() const;

  std::string name() const;

 protected:
  // Any lattice element that is negative gets directed to the
  // outer_universe. If outer_universe_index is also negative,
  // the particle is considered to be lost.
  std::vector<int32_t> lattice_universes;
  int32_t outer_universe_index;
  uint32_t id_;
  std::string name_;

};  // Lattice

#endif
