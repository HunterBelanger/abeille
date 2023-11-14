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
#ifndef LATTICE_UNIVERSE_H
#define LATTICE_UNIVERSE_H

#include <cstdint>
#include <geometry/lattice.hpp>
#include <geometry/universe.hpp>

#include <yaml-cpp/yaml.h>

#include <map>

//===========================================================================
// Externals from geometry
namespace geometry {
extern std::vector<std::shared_ptr<Universe>> universes;
extern std::vector<std::shared_ptr<Lattice>> lattices;
}  // namespace geometry

//===========================================================================
// Externsal from parser
extern std::map<uint32_t, size_t> universe_id_to_indx;
extern std::map<uint32_t, size_t> lattice_id_to_indx;

class LatticeUniverse : public Universe {
 public:
  LatticeUniverse(uint32_t i_lat, uint32_t i_id, std::string i_name);
  ~LatticeUniverse() = default;

  Cell* get_cell(Position r, Direction u, int32_t on_surf) const override;

  Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                 int32_t on_surf) const override;

  Boundary get_boundary_condition(const Position& r, const Direction& u,
                                  int32_t on_surf) const override final;

  Boundary lost_get_boundary(const Position& r, const Direction& u,
                             int32_t on_surf) const override final;

  bool contains_universe(uint32_t id) const override;

  uint64_t number_of_cell_instances(uint32_t id) const override;
  
  std::set<uint32_t> get_all_contained_cells() const override;

 private:
  uint32_t lattice_index;

};  // LatticeUniverse

//===========================================================================
// Non-Member functions
void make_lattice_universe(const YAML::Node& uni_node, const YAML::Node& input);

#endif
