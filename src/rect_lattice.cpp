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
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <geometry/surfaces/xplane.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <cmath>

RectLattice::RectLattice(uint32_t nx, uint32_t ny, uint32_t nz, double px,
                         double py, double pz, double xl, double yl, double zl,
                         uint32_t i_id, std::string i_name)
    : Lattice{i_id, i_name},
      Nx{nx},
      Ny{ny},
      Nz{nz},
      Px{px},
      Py{py},
      Pz{pz},
      Px_inv{1. / px},
      Py_inv{1. / py},
      Pz_inv{1. / pz},
      Xl{xl},
      Yl{yl},
      Zl{zl} {
  Xl = Xl - static_cast<double>(Nx) * 0.5 * Px;
  Yl = Yl - static_cast<double>(Ny) * 0.5 * Py;
  Zl = Zl - static_cast<double>(Nz) * 0.5 * Pz;
}

bool RectLattice::is_inside(Position r, Direction u) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  if ((nx < 0 || nx >= static_cast<int>(Nx)) ||
      (ny < 0 || ny >= static_cast<int>(Ny)) ||
      (nz < 0 || nz >= static_cast<int>(Nz))) {
    // Index is outside of lattice
    return false;
  } else {
    if (lattice_universes[linear_index(static_cast<uint32_t>(nx),
                                       static_cast<uint32_t>(ny),
                                       static_cast<uint32_t>(nz))] >= 0) {
      return true;
    } else {
      return false;
    }
  }
}

UniqueCell RectLattice::get_cell(Position r, Direction u,
                                 int32_t on_surf) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  UniqueCell ucell;

  if ((nx < 0 || nx >= static_cast<int>(Nx)) ||
      (ny < 0 || ny >= static_cast<int>(Ny)) ||
      (nz < 0 || nz >= static_cast<int>(Nz))) {
    // Index is outside of lattice, if outside_universe, try outer_universe
    if (outer_universe_index >= 0) {
      ucell =
          geometry::universes[static_cast<std::size_t>(outer_universe_index)]
              ->get_cell(r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    } else {
      // Location can not be found, return nullptr
      return ucell;
    }
  } else {
    const auto lin_indx =
        linear_index(static_cast<uint32_t>(nx), static_cast<uint32_t>(ny),
                     static_cast<uint32_t>(nz));
    if (lattice_universes[lin_indx] >= 0) {
      // Element is a valid fill, get cell from that universe
      // Transform coordinates to lattice elements locale frame
      Position r_local = r - tile_center(nx, ny, nz);
      const int32_t univ_indx = lattice_universes[lin_indx];
      ucell =
          geometry::universes[static_cast<std::size_t>(univ_indx)]->get_cell(
              r_local, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map[lin_indx].at(ucell.id);
      return ucell;
    } else {
      // Element is a dummy, try outer_universe
      if (outer_universe_index >= 0) {
        // outer_universe is give, get cell from that
        ucell =
            geometry::universes[static_cast<std::size_t>(outer_universe_index)]
                ->get_cell(r, u, on_surf);
        if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
        return ucell;
      } else {
        // No outer_universe provided, return nullptr
        return ucell;
      }
    }
  }
}

UniqueCell RectLattice::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                                 Direction u, int32_t on_surf) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  UniqueCell ucell;

  if ((nx < 0 || nx >= static_cast<int>(Nx)) ||
      (ny < 0 || ny >= static_cast<int>(Ny)) ||
      (nz < 0 || nz >= static_cast<int>(Nz))) {
    // Index is outside of lattice, if outside_universe, try outer_universe
    if (outer_universe_index >= 0) {
      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, true});

      // Go to outside universe
      ucell =
          geometry::universes[static_cast<std::size_t>(outer_universe_index)]
              ->get_cell(stack, r, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
      return ucell;
    } else {
      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

      // Location can not be found, return nullptr
      return ucell;
    }
  } else {
    const auto lin_indx =
        linear_index(static_cast<uint32_t>(nx), static_cast<uint32_t>(ny),
                     static_cast<uint32_t>(nz));
    if (lattice_universes[lin_indx] >= 0) {
      // Element is a valid fill, get cell from that universe
      // Transform coordinates to lattice elements locale frame
      Position r_local = r - tile_center(nx, ny, nz);
      const int32_t univ_indx = lattice_universes[lin_indx];

      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

      ucell =
          geometry::universes[static_cast<std::size_t>(univ_indx)]->get_cell(
              stack, r_local, u, on_surf);
      if (ucell) ucell.instance += cell_offset_map[lin_indx].at(ucell.id);
      return ucell;
    } else {
      // Element is a dummy, try outer_universe
      if (outer_universe_index >= 0) {
        // Save lattice info to stack
        stack.push_back(
            {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, true});

        // outer_universe is give, get cell from that
        ucell =
            geometry::universes[static_cast<std::size_t>(outer_universe_index)]
                ->get_cell(stack, r, u, on_surf);
        if (ucell) ucell.instance += cell_offset_map.back().at(ucell.id);
        return ucell;
      } else {
        // Save lattice info to stack
        stack.push_back(
            {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

        // No outer_universe provided, return nullptr
        return ucell;
      }
    }
  }
}

std::array<int32_t, 3> RectLattice::get_tile(Position r, Direction u) const {
  int32_t nx = static_cast<int32_t>(std::floor((r.x() - Xl) * Px_inv));
  int32_t ny = static_cast<int32_t>(std::floor((r.y() - Yl) * Py_inv));
  int32_t nz = static_cast<int32_t>(std::floor((r.z() - Zl) * Pz_inv));

  Position r_tile = tile_center(nx, ny, nz);

  // Bounds of tile. Must check all to see if we are on them, then use
  // direction to see if we are in or out of tile.
  double xl = r_tile.x() - Px * 0.5;
  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT && u.x() < 0.) nx--;

  double xh = r_tile.x() + Px * 0.5;
  if (std::abs(xh - r.x()) < SURFACE_COINCIDENT && u.x() >= 0.) nx++;

  double yl = r_tile.y() - Py * 0.5;
  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT && u.y() < 0.) ny--;

  double yh = r_tile.y() + Py * 0.5;
  if (std::abs(yh - r.y()) < SURFACE_COINCIDENT && u.y() >= 0.) ny++;

  double zl = r_tile.z() - Pz * 0.5;
  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT && u.z() < 0.) nz--;

  double zh = r_tile.z() + Pz * 0.5;
  if (std::abs(zh - r.z()) < SURFACE_COINCIDENT && u.z() >= 0.) nz++;

  return {nx, ny, nz};
}

double RectLattice::distance_to_tile_boundary(
    Position r_local, Direction u, std::array<int32_t, 3> tile) const {
  // Get the position of the center of the tile
  Position center = tile_center(tile[0], tile[1], tile[2]);

  // Position relatice to center of tile
  Position r_tile = r_local - center;

  double dist = INF;

  // Check all surfaces
  const double diff_xl = -Px * 0.5 - r_tile.x();
  const double diff_xh = Px * 0.5 - r_tile.x();
  const double diff_yl = -Py * 0.5 - r_tile.y();
  const double diff_yh = Py * 0.5 - r_tile.y();
  const double diff_zl = -Pz * 0.5 - r_tile.z();
  const double diff_zh = Pz * 0.5 - r_tile.z();

  const double ux_inv = 1. / u.x();
  const double uy_inv = 1. / u.y();
  const double uz_inv = 1. / u.z();

  const double d_xl = diff_xl * ux_inv;
  const double d_xh = diff_xh * ux_inv;
  const double d_yl = diff_yl * uy_inv;
  const double d_yh = diff_yh * uy_inv;
  const double d_zl = diff_zl * uz_inv;
  const double d_zh = diff_zh * uz_inv;

  if (d_xl > 0. && d_xl < dist && std::abs(diff_xl) > 100 * SURFACE_COINCIDENT)
    dist = d_xl;
  if (d_xh > 0. && d_xh < dist && std::abs(diff_xh) > 100 * SURFACE_COINCIDENT)
    dist = d_xh;
  if (d_yl > 0. && d_yl < dist && std::abs(diff_yl) > 100 * SURFACE_COINCIDENT)
    dist = d_yl;
  if (d_yh > 0. && d_yh < dist && std::abs(diff_yh) > 100 * SURFACE_COINCIDENT)
    dist = d_yh;
  if (d_zl > 0. && d_zl < dist && std::abs(diff_zl) > 100 * SURFACE_COINCIDENT)
    dist = d_zl;
  if (d_zh > 0. && d_zh < dist && std::abs(diff_zh) > 100 * SURFACE_COINCIDENT)
    dist = d_zh;

  return dist;
}

void RectLattice::set_elements(std::vector<int32_t> univs) {
  if (univs.size() == (Nx * Ny * Nz)) {
    lattice_universes = univs;
  } else {
    // Wrong number of elements provided, throw error
    fatal_error("Improper number of elements provided to lattice.");
  }
}

size_t RectLattice::linear_index(uint32_t nx, uint32_t ny, uint32_t nz) const {
  if (nz >= Nz || nx >= Nx || ny >= Ny) {
    // Bad index
    fatal_error("Invalid lattice indecies");
  }

  size_t indx;
  indx = static_cast<size_t>(nz * (Nx * Ny) + nx * Ny + ny);
  return indx;
}

Position RectLattice::tile_center(int nx, int ny, int nz) const {
  double x = (static_cast<double>(nx) + 0.5) * Px + Xl;
  double y = (static_cast<double>(ny) + 0.5) * Py + Yl;
  double z = (static_cast<double>(nz) + 0.5) * Pz + Zl;

  return Position(x, y, z);
}

void make_rect_lattice(const YAML::Node& latt_node, const YAML::Node& input) {
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
  shape.resize(3);
  if (latt_node["shape"] && latt_node["shape"].IsSequence() &&
      latt_node["shape"].size() == 3) {
    for (size_t s = 0; s < 3; s++) {
      shape[s] = latt_node["shape"][s].as<uint32_t>();
    }
  } else {
    fatal_error("Lattice must have a valid shape.");
  }

  // Get pitch
  std::vector<double> pitch;
  pitch.resize(3);
  if (latt_node["pitch"] && latt_node["pitch"].IsSequence() &&
      latt_node["pitch"].size() == 3) {
    for (size_t s = 0; s < 3; s++) {
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
    if (latt_node["universes"].size() == shape[0] * shape[1] * shape[2]) {
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
      fatal_error("Lattice instance has improper number of universes.");
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

  // Make lattice
  std::shared_ptr<Lattice> lat = std::make_shared<RectLattice>(
      shape[0], shape[1], shape[2], pitch[0], pitch[1], pitch[2], origin[0],
      origin[1], origin[2], id, name);

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
