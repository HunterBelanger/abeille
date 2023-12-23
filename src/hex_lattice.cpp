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
#include <geometry/hex_lattice.hpp>
#include <utils/error.hpp>

#include <cmath>

HexLattice::HexLattice(uint32_t nrings, uint32_t nz, double p, double pz,
                       double x, double y, double z, Top t, uint32_t i_id,
                       std::string i_name)
    : Lattice{i_id, i_name},
      pitch_{p},
      pitch_z_{pz},
      X_o{x},
      Y_o{y},
      Z_o{z},
      Nrings{nrings},
      Nz{nz},
      Nhex{0},
      width{0},
      mid_qr{0},
      top_{t} {
  Nhex = 0;
  for (uint32_t r = 0; r < Nrings; r++) {
    if (r == 0)
      Nhex += 1;
    else
      Nhex += 6 * r;
  }

  width = 2 * (Nrings - 1) + 1;

  mid_qr = width / 2;
}

bool HexLattice::is_inside(Position r, Direction /*u*/) const {
  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 2> qr = get_nearest_hex(r_o);

  // Get the ring of hex tile
  uint32_t ring = get_ring(qr);

  // See if valid ring or not
  if (ring >= Nrings)
    // Invalid ring
    return false;

  // Check z bin now
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  int32_t nz = static_cast<int32_t>(std::floor((r.z() - Z_low) / pitch_z_));
  if (nz < 0 || nz >= static_cast<int32_t>(Nz)) return false;

  return true;
}

UniqueCell HexLattice::get_cell(Position r, Direction u,
                                int32_t on_surf) const {
  UniqueCell ucell;

  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 3> qrz = get_tile(r_o, u);

  // Get the ring of hex tile
  uint32_t ring = get_ring({qrz[0], qrz[1]});

  // See if valid ring or not
  if (ring >= Nrings) {
    // Invalid ring
    if (outer_universe_index == -1) {
      return ucell;
    } else {
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Check z bin now
  if (qrz[2] < 0 || qrz[2] >= static_cast<int32_t>(Nz)) {
    if (outer_universe_index == -1) {
      return ucell;
    } else {
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Inside a lattice bin
  const size_t indx = linear_index({qrz[0], qrz[1]}, qrz[2]);

  // If -1, send to outside universe
  if (lattice_universes[indx] == -1) {
    if (outer_universe_index == -1) {
      return ucell;
    } else {
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Move coordinates to tile center
  Position center = tile_center(qrz[0], qrz[1], qrz[2]);
  Position r_tile = r_o - center;
  ucell = geometry::universes[static_cast<uint32_t>(lattice_universes[indx])]
              ->get_cell(r_tile, u, on_surf);
  if (ucell) ucell.instance += cell_offset_map[indx].at(ucell.id);
  return ucell;
}

UniqueCell HexLattice::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                                Direction u, int32_t on_surf) const {
  UniqueCell ucell;

  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 3> qrz = get_tile(r_o, u);

  // Get the ring of hex tile
  uint32_t ring = get_ring({qrz[0], qrz[1]});

  // See if valid ring or not
  if (ring >= Nrings) {
    // Invalid ring
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return ucell;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(stack, r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Check z bin now
  if (qrz[2] < 0 || qrz[2] >= static_cast<int32_t>(Nz)) {
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return ucell;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(stack, r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Inside a lattice bin
  const size_t indx = linear_index({qrz[0], qrz[1]}, qrz[2]);

  // If -1, send to outside universe
  if (lattice_universes[indx] == -1) {
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return ucell;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      ucell = geometry::universes[static_cast<uint32_t>(outer_universe_index)]
                  ->get_cell(stack, r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    }
  }

  // Move coordinates to tile center
  Position center = tile_center(qrz[0], qrz[1], qrz[2]);
  Position r_tile = r_o - center;
  // Save info to stack
  stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
  ucell = geometry::universes[static_cast<uint32_t>(lattice_universes[indx])]
              ->get_cell(stack, r_tile, u, on_surf);
  if (ucell) ucell.instance += cell_offset_map[indx].at(ucell.id);
  return ucell;
}

void HexLattice::set_elements(std::vector<int32_t> univs) {
  // Make sure proper number of universes are provided
  if (univs.size() != Nhex * Nz) {
    fatal_error("Improper number of universes for HexLattice.");
  }
  lattice_universes.resize(width * width * Nz, -1);

  size_t indx = 0;
  for (uint32_t az = 0; az < Nz; az++) {
    for (uint32_t ar = 0; ar < width; ar++) {
      for (uint32_t aq = 0; aq < width; aq++) {
        int32_t q = static_cast<int32_t>(aq - mid_qr);
        int32_t r = static_cast<int32_t>(ar - mid_qr);

        uint32_t ring = get_ring({q, r});

        if (ring < Nrings) {
          lattice_universes[linear_index({q, r}, static_cast<int32_t>(az))] =
              univs[indx];
          indx += 1;
        }
      }
    }
  }
}

size_t HexLattice::linear_index(std::array<int32_t, 2> qr, int32_t z) const {
  // In qr coordinates shifted by mid_qr, both q and r can range from 0 to
  // width - 1. This then forms a width*width square array.
  size_t indx = 0;
  uint32_t q = static_cast<uint32_t>(qr[0]) + mid_qr;
  uint32_t r = static_cast<uint32_t>(qr[1]) + mid_qr;
  indx = static_cast<uint32_t>(z) * (width * width) + r * (width) + q;

  if (indx >= lattice_universes.size()) {
    fatal_error("Invalid index for HexLattice.");
  }

  return indx;
}

std::array<int32_t, 2> HexLattice::get_nearest_hex(Position p) const {
  double q, r, det;
  if (top_ == Top::Flat) {
    det = -pitch_ * pitch_ * sin_pi_3;

    q = (pitch_ / det) * (0. * p.x() - 1. * p.y());
    r = (pitch_ / det) * (-sin_pi_3 * p.x() + cos_pi_3 * p.y());
  } else {
    det = -pitch_ * pitch_ * cos_pi_6;

    q = (pitch_ / det) * (sin_pi_6 * p.x() - cos_pi_6 * p.y());
    r = (pitch_ / det) * (-1. * p.x() - 0. * p.y());
  }

  // Do rounding magic to get nearest hexagon
  // First convert to cube coordinates
  double x = q;
  double z = r;
  double y = -x - z;

  // Do rounding
  double rx = std::round(x);
  double ry = std::round(y);
  double rz = std::round(z);

  // Differences
  double x_diff = std::abs(rx - x);
  double y_diff = std::abs(ry - y);
  double z_diff = std::abs(rz - z);

  // Check three cases, correct rounded coordinates
  if (x_diff > y_diff && x_diff > z_diff) {
    rx = -ry - rz;
  } else if (y_diff > x_diff && y_diff > z_diff) {
    ry = -rx - rz;
  } else {
    rz = -rx - ry;
  }

  // Convert back to axial coordinates for return
  return {static_cast<int32_t>(rx), static_cast<int32_t>(rz)};
}

std::array<int32_t, 3> HexLattice::get_tile(Position p, Direction /*u*/) const {
  std::array<int32_t, 2> qr = get_nearest_hex(p);
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  int32_t nz = static_cast<int32_t>(std::floor((p.z() - Z_low) / pitch_z_));
  return {qr[0], qr[1], nz};
}

Position HexLattice::get_hex_center(std::array<int32_t, 2> qr) const {
  double x, y;
  double q = static_cast<double>(qr[0]);
  double r = static_cast<double>(qr[1]);

  if (top_ == Top::Flat) {
    x = pitch_ * (cos_pi_3 * q + 1. * r);
    y = pitch_ * (sin_pi_3 * q + 0. * r);
  } else {
    x = pitch_ * (0. * q + cos_pi_6 * r);
    y = pitch_ * (1. * q + sin_pi_6 * r);
  }

  return Position(x, y, 0.);
}

Position HexLattice::tile_center(int q, int r, int nz) const {
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  Position r_hex = get_hex_center({q, r});
  double z = (static_cast<double>(nz) + 0.5) * pitch_z_ + Z_low;
  return Position(r_hex.x(), r_hex.y(), z);
}

uint32_t HexLattice::get_ring(std::array<int32_t, 2> qr) const {
  int32_t x = qr[0];
  int32_t z = qr[1];
  int32_t y = -x - z;

  uint32_t ax = static_cast<uint32_t>(std::abs(x));
  uint32_t ay = static_cast<uint32_t>(std::abs(y));
  uint32_t az = static_cast<uint32_t>(std::abs(z));

  // Ring is the max abs of x, y z
  uint32_t temp_abs = std::max(ax, ay);
  return std::max(temp_abs, az);
}

double HexLattice::distance_to_tile_boundary(
    Position r_local, Direction u, std::array<int32_t, 3> tile) const {
  Position center = tile_center(tile[0], tile[1], tile[2]);
  Position r_tile = r_local - center;

  if (top_ == Top::Pointy)
    return distance_to_tile_boundary_pointy(r_tile, u);
  else
    return distance_to_tile_boundary_flat(r_tile, u);
}

double HexLattice::distance_to_tile_boundary_pointy(Position r_tile,
                                                    Direction u) const {
  // Points needed to calculate slopes for lines. Only need two sets
  // because of symmetry reasons.
  //    (x1,y1)
  /*      /\
         |  |(x2,y2)  */
  double x1 = 0.;
  double y1 = pitch_ / (2. * cos_pi_6);
  double x2 = pitch_ / 2.;
  double y2 = pitch_ / (2. * sin_pi_6);

  double d = INF;

  // Get distance to all six hexagon surfaces
  double d1 = distance_to_line(r_tile, u, x1, y1, x2, y1);
  double d2 = distance_to_line(r_tile, u, x2, y2, x2, -y2);
  double d3 = distance_to_line(r_tile, u, x2, -y2, x1, -y1);
  double d4 = distance_to_line(r_tile, u, x1, -y1, -x2, -y2);
  double d5 = distance_to_line(r_tile, u, -x2, -y2, -x2, y2);
  double d6 = distance_to_line(r_tile, u, -x2, y2, x1, y1);

  // Get distance to z planes
  double dzl = (-pitch_z_ * 0.5 - r_tile.z()) / u.z();
  double dzu = (pitch_z_ * 0.5 - r_tile.z()) / u.z();

  if (d1 > 0. && d1 < d) d = d1;
  if (d2 > 0. && d2 < d) d = d2;
  if (d3 > 0. && d3 < d) d = d3;
  if (d4 > 0. && d4 < d) d = d4;
  if (d5 > 0. && d5 < d) d = d5;
  if (d6 > 0. && d6 < d) d = d6;
  if (dzl > 0. && dzl < d) d = dzl;
  if (dzu > 0. && dzu < d) d = dzu;

  return d;
}

double HexLattice::distance_to_tile_boundary_flat(Position r_tile,
                                                  Direction u) const {
  // Points needed to calculate slopes for lines. Only need two sets
  // because of symmetry reasons.
  //
  /*        --- (x1, y1)
           /   \
                 (x2, y2)  */
  double x1 = pitch_ / (2. * sin_pi_6);
  double y1 = pitch_ / 2.;
  double x2 = pitch_ / (2. * cos_pi_6);
  double y2 = 0.;

  double d = INF;

  // Get distance to all six hexagon surfaces
  double d1 = distance_to_line(r_tile, u, x1, y1, x2, y1);
  double d2 = distance_to_line(r_tile, u, x2, y2, x1, -y1);
  double d3 = distance_to_line(r_tile, u, x1, -y1, -x1, -y1);
  double d4 = distance_to_line(r_tile, u, -x1, -y1, -x2, -y2);
  double d5 = distance_to_line(r_tile, u, -x2, -y2, -x1, y1);
  double d6 = distance_to_line(r_tile, u, -x1, y1, x1, y1);

  // Get distance to z planes
  double dzl = (-pitch_z_ * 0.5 - r_tile.z()) / u.z();
  double dzu = (pitch_z_ * 0.5 - r_tile.z()) / u.z();

  if (d1 > 0. && d1 < d) d = d1;
  if (d2 > 0. && d2 < d) d = d2;
  if (d3 > 0. && d3 < d) d = d3;
  if (d4 > 0. && d4 < d) d = d4;
  if (d5 > 0. && d5 < d) d = d5;
  if (d6 > 0. && d6 < d) d = d6;
  if (dzl > 0. && dzl < d) d = dzl;
  if (dzu > 0. && dzu < d) d = dzu;

  return d;
}

double HexLattice::distance_to_line(Position r, Direction u, double x1,
                                    double y1, double x2, double y2) const {
  // Get A, B, and D for plane from y - y1 = m (x - x1) ; m = (y2 - y1)/(x2 -
  // x1)
  double A = y2 - y1;
  double B = x1 - x2;
  double D = (x2 - x1) * y1 - (y2 - y1) * x1;

  double num = D - A * r.x() - B * r.y();
  double denom = A * u.x() + B * u.y();

  double d = num / denom;

  if (d < 0.)
    return INFINITY;
  else
    return d;
}

void make_hex_lattice(const YAML::Node& latt_node, const YAML::Node& input) {
  // Get id
  uint32_t id = 0;
  if (latt_node["id"] && latt_node["id"].IsScalar()) {
    id = latt_node["id"].as<uint32_t>();
  } else {
    fatal_error("Lattice must have a valid id.");
  }

  // Get name if present
  std::string name = "";
  if (latt_node["name"] && latt_node["name"].IsScalar()) {
    name = latt_node["name"].as<std::string>();
  }

  // Get shape
  std::vector<uint32_t> shape;
  shape.resize(2);
  if (latt_node["shape"] && latt_node["shape"].IsSequence() &&
      latt_node["shape"].size() == 2) {
    for (size_t s = 0; s < 2; s++) {
      shape[s] = latt_node["shape"][s].as<uint32_t>();
    }
  } else {
    fatal_error("Lattice must have a valid shape.");
  }

  // Get pitch
  std::vector<double> pitch;
  pitch.resize(2);
  if (latt_node["pitch"] && latt_node["pitch"].IsSequence() &&
      latt_node["pitch"].size() == 2) {
    for (size_t s = 0; s < 2; s++) {
      pitch[s] = latt_node["pitch"][s].as<double>();
    }
  } else {
    fatal_error("Lattice must have a valid pitch.");
  }

  // Get origin
  std::vector<double> origin;
  origin.resize(3);
  if (latt_node["origin"] && latt_node["origin"].IsSequence() &&
      latt_node["origin"].size() == 3) {
    for (size_t s = 0; s < 3; s++) {
      origin[s] = latt_node["origin"][s].as<double>();
    }
  } else {
    fatal_error("Lattice must have a valid origin.");
  }

  // Vector for universes
  std::vector<int32_t> uni_indicies;
  if (latt_node["universes"] && latt_node["universes"].IsSequence()) {
    // Go through and check each universe
    for (size_t u = 0; u < latt_node["universes"].size(); u++) {
      int32_t u_id = latt_node["universes"][u].as<int32_t>();
      if (u_id == -1) {
        uni_indicies.push_back(u_id);
      } else if (universe_id_to_indx.find(static_cast<uint32_t>(u_id)) ==
                 universe_id_to_indx.end()) {
        // Need to find universe
        find_universe(input, static_cast<uint32_t>(u_id));
        uni_indicies.push_back(static_cast<int32_t>(
            universe_id_to_indx[static_cast<uint32_t>(u_id)]));
      } else {
        uni_indicies.push_back(static_cast<int32_t>(
            universe_id_to_indx[static_cast<uint32_t>(u_id)]));
      }
    }
  } else {
    fatal_error("Lattice instance must have a valid universes definition.");
  }

  // Make sure id not taken
  if (universe_id_to_indx.find(id) != universe_id_to_indx.end()) {
    std::stringstream mssg;
    mssg << "Universe id " << id << " appears multiple times.";
    fatal_error(mssg.str());
  }

  // Get top type
  HexLattice::Top top = HexLattice::Top::Pointy;
  if (latt_node["top"] && latt_node["top"].IsScalar()) {
    if (latt_node["top"].as<std::string>() == "pointy") {
      top = HexLattice::Top::Pointy;
    } else if (latt_node["top"].as<std::string>() == "flat") {
      top = HexLattice::Top::Flat;
    } else {
      std::stringstream mssg;
      mssg << " Uknown top for hexagonal lattice "
           << latt_node["top"].as<std::string>() << ".";
      fatal_error(mssg.str());
    }
  }

  // Make lattice
  std::shared_ptr<Lattice> lat = std::make_shared<HexLattice>(
      shape[0], shape[1], pitch[0], pitch[1], origin[0], origin[1], origin[2],
      top, id, name);

  // Set universes
  lat->set_elements(uni_indicies);

  // Get outside
  int32_t out_id;
  if (latt_node["outer"] && latt_node["outer"].IsScalar()) {
    out_id = latt_node["outer"].as<int32_t>();
    if (out_id != -1) {
      // Find outside universe
      if (universe_id_to_indx.find(static_cast<uint32_t>(out_id)) ==
          universe_id_to_indx.end()) {
        // Need to find universe
        find_universe(input, static_cast<uint32_t>(out_id));
      }
      lat->set_outisde_universe(static_cast<int32_t>(
          universe_id_to_indx[static_cast<uint32_t>(out_id)]));
    }
  }

  // Set lattice
  universe_id_to_indx[id] = geometry::universes.size();
  geometry::universes.push_back(lat);
}
