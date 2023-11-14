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
#ifndef CELL_UNIVERSE_H
#define CELL_UNIVERSE_H

#include <geometry/universe.hpp>

#include <yaml-cpp/yaml.h>

#include <map>

//===========================================================================
// Externals from geometry
namespace geometry {
extern std::vector<std::shared_ptr<Cell>> cells;
extern std::vector<std::shared_ptr<Universe>> universes;
}  // namespace geometry

//===========================================================================
// Externals from parser
extern std::map<uint32_t, size_t> cell_id_to_indx;
extern std::map<uint32_t, size_t> universe_id_to_indx;

//===========================================================================
// CellUniverse class
class CellUniverse : public Universe {
 public:
  CellUniverse(std::vector<uint32_t> i_ind, uint32_t i_id, std::string i_name);
  ~CellUniverse() = default;

  Cell* get_cell(Position r, Direction u, int32_t on_surf) const override;

  Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                 int32_t on_surf) const override;

  Boundary get_boundary_condition(const Position& r, const Direction& u,
                                  int32_t on_surf) const override final;

  Boundary lost_get_boundary(const Position& r, const Direction& u,
                             int32_t on_surf) const override final;

  bool contains_universe(uint32_t id) const override;

 private:
  std::vector<uint32_t> cell_indicies;

};  // CellUniverse

//===========================================================================
// Non-Member functions
void make_cell_universe(const YAML::Node& uni_node);

#endif
