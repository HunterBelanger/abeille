/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#include <cmath>
#include <geometry/rect_lattice.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

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
    if (lattice_universes[linear_index(nx, ny, nz)] >= 0) {
      return true;
    } else {
      return false;
    }
  }
}

Cell* RectLattice::get_cell(Position r, Direction u, int32_t on_surf) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  if ((nx < 0 || nx >= static_cast<int>(Nx)) ||
      (ny < 0 || ny >= static_cast<int>(Ny)) ||
      (nz < 0 || nz >= static_cast<int>(Nz))) {
    // Index is outside of lattice, if outside_universe, try outer_universe
    if (outer_universe_index >= 0) {
      return geometry::universes[outer_universe_index]->get_cell(r, u, on_surf);
    } else {
      // Location can not be found, return nullptr
      return nullptr;
    }
  } else {
    if (lattice_universes[linear_index(nx, ny, nz)] >= 0) {
      // Element is a valid fill, get cell from that universe
      // Transform coordinates to lattice elements locale frame
      Position r_local = r - tile_center(nx, ny, nz);
      int32_t univ_indx = lattice_universes[linear_index(nx, ny, nz)];
      return geometry::universes[univ_indx]->get_cell(r_local, u, on_surf);
    } else {
      // Element is a dummy, try outer_universe
      if (outer_universe_index >= 0) {
        // outer_universe is give, get cell from that
        return geometry::universes[outer_universe_index]->get_cell(r, u,
                                                                   on_surf);
      } else {
        // No outer_universe provided, return nullptr
        return nullptr;
      }
    }
  }
}

Cell* RectLattice::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                            Direction u, int32_t on_surf) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  if ((nx < 0 || nx >= static_cast<int>(Nx)) ||
      (ny < 0 || ny >= static_cast<int>(Ny)) ||
      (nz < 0 || nz >= static_cast<int>(Nz))) {
    // Index is outside of lattice, if outside_universe, try outer_universe
    if (outer_universe_index >= 0) {
      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, true});

      // Go to outside universe
      return geometry::universes[outer_universe_index]->get_cell(stack, r, u,
                                                                 on_surf);
    } else {
      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

      // Location can not be found, return nullptr
      return nullptr;
    }
  } else {
    if (lattice_universes[linear_index(nx, ny, nz)] >= 0) {
      // Element is a valid fill, get cell from that universe
      // Transform coordinates to lattice elements locale frame
      Position r_local = r - tile_center(nx, ny, nz);
      int32_t univ_indx = lattice_universes[linear_index(nx, ny, nz)];

      // Save lattice info to stack
      stack.push_back(
          {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

      return geometry::universes[univ_indx]->get_cell(stack, r_local, u,
                                                      on_surf);
    } else {
      // Element is a dummy, try outer_universe
      if (outer_universe_index >= 0) {
        // Save lattice info to stack
        stack.push_back(
            {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, true});

        // outer_universe is give, get cell from that
        return geometry::universes[outer_universe_index]->get_cell(stack, r, u,
                                                                   on_surf);
      } else {
        // Save lattice info to stack
        stack.push_back(
            {GeoLilyPad::PadType::Lattice, id_, r, {nx, ny, nz}, false});

        // No outer_universe provided, return nullptr
        return nullptr;
      }
    }
  }
}

std::array<int32_t, 3> RectLattice::get_tile(Position r, Direction u) const {
  int32_t nx = static_cast<int32_t>(std::floor((r.x() - Xl) / Px));
  int32_t ny = static_cast<int32_t>(std::floor((r.y() - Yl) / Py));
  int32_t nz = static_cast<int32_t>(std::floor((r.z() - Zl) / Pz));

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
  // Check if we are inside the lattice, or outside. If we are outside,
  // we need to do something different.
  if (tile[0] < 0 || tile[1] < 0 || tile[2] < 0 ||
      tile[0] >= static_cast<int32_t>(Nx) ||
      tile[1] >= static_cast<int32_t>(Ny) ||
      tile[2] > static_cast<int32_t>(Nz)) {
    double dmin = (Xl - r_local.x()) / u.x();
    double dmax = (Xl+(Nx*Px) - r_local.x()) / u.x();
    if (dmin > dmax) std::swap(dmin, dmax);

    double dymin = (Yl - r_local.y()) / u.y();
    double dymax = (Yl+(Ny*Py) - r_local.y()) / u.y();
    if (dymin > dymax) std::swap(dymin, dymax);

    if ((dmin > dymax) || (dymin > dmax)) return INF;

    if (dymin > dmin) dmin = dymin;

    if (dymax < dmax) dmax = dymax; 

    double dzmin = (Zl - r_local.z()) / u.z();
    double dzmax = (Zl+(Nz*Pz) - r_local.z()) / u.z();
    if (dzmin > dzmax) std::swap(dzmin, dzmax);

    if ((dmin > dzmax) || (dzmin > dmax)) return INF;

    if (dzmin > dmin) dmin = dzmin;

    if (dzmax < dmax) dmax = dzmax;

    if (dmin < 0. && dmax < 0.) return INF;

    if (dmin < 0.) return dmax;

    return dmin;
  }

  // If we are here, we are actually inside a lattice tile,
  // so we just find the distance to the next tile.

  // Get the position of the center of the tile
  Position center = tile_center(tile[0], tile[1], tile[2]);

  // Position relatice to center of tile
  Position r_tile = r_local - center;

  double dist = INF;

  // Check all surfaces
  double diff_xl = -Px * 0.5 - r_tile.x();
  double diff_xh = Px * 0.5 - r_tile.x();
  double diff_yl = -Py * 0.5 - r_tile.y();
  double diff_yh = Py * 0.5 - r_tile.y();
  double diff_zl = -Pz * 0.5 - r_tile.z();
  double diff_zh = Pz * 0.5 - r_tile.z();

  double d_xl = diff_xl / u.x();
  double d_xh = diff_xh / u.x();
  double d_yl = diff_yl / u.y();
  double d_yh = diff_yh / u.y();
  double d_zl = diff_zl / u.z();
  double d_zh = diff_zh / u.z();
  
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
    std::string mssg = "Improper number of elements provided to lattice.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

size_t RectLattice::linear_index(uint32_t nx, uint32_t ny, uint32_t nz) const {
  if (nz >= Nz || nx >= Nx || ny >= Ny) {
    // Bad index
    std::string mssg = "Invalid lattice indecies";
    fatal_error(mssg, __FILE__, __LINE__);
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

void make_rect_lattice(YAML::Node latt_node, YAML::Node input) {
  // Get id
  uint32_t id = 0;
  if (latt_node["id"] && latt_node["id"].IsScalar()) {
    id = latt_node["id"].as<uint32_t>();
  } else {
    std::string mssg = "Lattice must have a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
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
    std::string mssg = "Lattice must have a valid shape.";
    fatal_error(mssg, __FILE__, __LINE__);
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
    std::string mssg = "Lattice must have a valid pitch.";
    fatal_error(mssg, __FILE__, __LINE__);
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
    std::string mssg = "Lattice must have a valid origin.";
    fatal_error(mssg, __FILE__, __LINE__);
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
        } else if (universe_id_to_indx.find(u_id) ==
                   universe_id_to_indx.end()) {
          // Need to find universe
          find_universe(input, u_id);
          uni_indicies.push_back(
              static_cast<int32_t>(universe_id_to_indx[u_id]));
        } else {
          uni_indicies.push_back(
              static_cast<int32_t>(universe_id_to_indx[u_id]));
        }
      }
    } else {
      std::string mssg = "Lattice instance has improper number of universes.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg =
        "Lattice instance must have a valid universes definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure lattice id not taken
  if (lattice_id_to_indx.find(id) != lattice_id_to_indx.end()) {
    std::string mssg = "Lattice id " + std::to_string(id) +
                       " appears"
                       " multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
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
      if (universe_id_to_indx.find(out_id) == universe_id_to_indx.end()) {
        // Need to find universe
        find_universe(input, out_id);
      }
      lat->set_outisde_universe(
          static_cast<int32_t>(universe_id_to_indx[out_id]));
    }
  }

  // Set lattice
  lattice_id_to_indx[id] = geometry::lattices.size();
  geometry::lattices.push_back(lat);
}
