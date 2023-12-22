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
#ifndef CELL_H
#define CELL_H

#include <geometry/surfaces/surface.hpp>
#include <materials/material.hpp>

#include <yaml-cpp/yaml.h>

#include <limits>
#include <map>
#include <memory>
#include <vector>

class Cell;

//============================================================================
// Operators for rpn cell definitions
enum OP : int32_t {
  L_PAR = std::numeric_limits<int32_t>::max(),
  R_PAR = std::numeric_limits<int32_t>::max() - 1,
  COMP = std::numeric_limits<int32_t>::max() - 2,
  INTR = std::numeric_limits<int32_t>::max() - 3,
  UNIN = std::numeric_limits<int32_t>::max() - 4
};

//============================================================================
// Externals from geometry.hpp
namespace geometry {
extern std::vector<std::shared_ptr<Surface>> surfaces;
extern std::vector<std::shared_ptr<Cell>> cells;
}  // namespace geometry

//============================================================================
// Universe class definition
class Universe;

//============================================================================
// Eternals from parser.hpp
extern std::map<uint32_t, size_t> surface_id_to_indx;
extern std::map<uint32_t, size_t> cell_id_to_indx;

//============================================================================
// Cell Class
class Cell {
 public:
  enum class Fill { Material, Universe };

 public:
  Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Material> material,
       uint32_t i_id, std::string i_name);
  Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Universe> universe,
       uint32_t i_id, std::string i_name);
  ~Cell() = default;

  Fill fill() const { return fill_; }

  bool vacuum_or_reflective() const { return vacuum_or_reflective_; }

  bool is_inside(const Position& r, const Direction& u, int32_t on_surf) const;

  std::pair<double, int32_t> distance_to_boundary(const Position& r,
                                                  const Direction& u,
                                                  int32_t on_surf) const;

  std::pair<double, int32_t> distance_to_boundary_condition(
      const Position& r, const Direction& u, int32_t on_surf) const;

  Material* material() { return material_raw_; }

  Universe* universe() { return universe_raw_; }

  uint32_t id() const;

  const std::string& name() const;

 private:
  Fill fill_;
  bool simple = true;
  bool vacuum_or_reflective_ = false;
  std::vector<int32_t> rpn;  // Surface definition of cell
  uint32_t id_;
  std::string name_;

  bool is_inside_simple(const Position& r, const Direction& u,
                        int32_t on_surf) const;
  bool is_inside_complex(const Position& r, const Direction& u,
                         int32_t on_surf) const;

  void check_for_bc();
  void simplify();

  std::shared_ptr<Material> material_;
  Material* material_raw_;

  std::shared_ptr<Universe> universe_;
  Universe* universe_raw_;
};  // Cell

//============================================================================
// UniqueCell
struct UniqueCell {
  Cell* cell = nullptr;
  uint32_t id = 0;
  uint32_t instance = 0;

  operator bool() const { return cell != nullptr; }
};

//===========================================================================
// Non-Member Functions
std::vector<int32_t> infix_to_rpn(const std::vector<int32_t>& infix);

void make_cell(const YAML::Node& cell_node, const YAML::Node& input);

#endif
