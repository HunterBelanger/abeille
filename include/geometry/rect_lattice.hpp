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
#ifndef RECT_LATTICE_H
#define RECT_LATTICE_H

#include <geometry/lattice.hpp>

class RectLattice : public Lattice {
 public:
  RectLattice(uint32_t nx, uint32_t ny, uint32_t nz, double px, double py,
              double pz, double xl, double yl, double zl, uint32_t i_id,
              std::string i_name);
  ~RectLattice() = default;

  bool is_inside(Position r, Direction u) const override;

  std::array<int32_t, 3> get_tile(Position r, Direction u) const override;

  UniqueCell get_cell(Position r, Direction u, int32_t on_surf) const override;

  UniqueCell get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                      int32_t on_surf) const override;

  double distance_to_tile_boundary(Position r_local, Direction u,
                                   std::array<int32_t, 3> tile) const override;

  void set_elements(std::vector<int32_t> univs) override;

 private:
  // Values required to describe a rectilinear lattice
  uint32_t Nx, Ny, Nz;            // Number of elements along each axis
  double Px, Py, Pz;              // Pitch along each axis
  double Px_inv, Py_inv, Pz_inv;  // Pitch along each axis
  double Xl, Yl, Zl;  // Lowest value of lattice along each coordinate

  size_t linear_index(uint32_t nx, uint32_t ny, uint32_t nz) const;

  Position tile_center(int nx, int ny, int nz) const;

};  // RectLattice

void make_rect_lattice(const YAML::Node& latt_node, const YAML::Node& input);

#endif
