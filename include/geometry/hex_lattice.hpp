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
#ifndef HEX_LATTICE_H
#define HEX_LATTICE_H

#include <geometry/lattice.hpp>
#include <utils/constants.hpp>

#include <yaml-cpp/yaml.h>

#include <array>

class HexLattice : public Lattice {
 public:
  enum Top { Pointy, Flat };

  HexLattice(uint32_t nrings, uint32_t nz, double p, double pz, double x,
             double y, double z, Top t, uint32_t i_id, std::string i_name);
  ~HexLattice() = default;

  bool is_inside(Position r, Direction u) const override;

  std::array<int32_t, 3> get_tile(Position p, Direction u) const override;

  UniqueCell get_cell(Position r, Direction u, int32_t on_surf) const override;

  UniqueCell get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                      int32_t on_surf) const override;

  double distance_to_tile_boundary(Position r_local, Direction u,
                                   std::array<int32_t, 3> tile) const override;

  void set_elements(std::vector<int32_t> univs) override;

 private:
  double pitch_, pitch_z_;
  double X_o, Y_o, Z_o;  // Origin of center hex
  const double cos_pi_6 = std::cos(PI / 6.0);
  const double sin_pi_6 = std::sin(PI / 6.0);
  const double cos_pi_3 = std::cos(PI / 3.0);
  const double sin_pi_3 = std::sin(PI / 3.0);
  uint32_t Nrings, Nz, Nhex;
  uint32_t width, mid_qr;
  Top top_;

  std::array<int32_t, 2> get_nearest_hex(Position p) const;

  Position get_hex_center(std::array<int32_t, 2> qr) const;

  Position tile_center(int q, int r, int nz) const;

  uint32_t get_ring(std::array<int32_t, 2> qr) const;

  double distance_to_tile_boundary_pointy(Position r_tile, Direction u) const;

  double distance_to_tile_boundary_flat(Position r_tile, Direction u) const;

  double distance_to_line(Position r, Direction u, double x1, double y1,
                          double x2, double y2) const;

  size_t linear_index(std::array<int32_t, 2> qr, int32_t z) const;

};  // HexLattice

void make_hex_lattice(const YAML::Node& latt_node, const YAML::Node& input);

#endif
