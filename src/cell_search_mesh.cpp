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
#include <geometry/cell_search_mesh.hpp>
#include <geometry/geometry.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

extern std::map<uint32_t, size_t> cell_id_to_indx;

CellSearchMesh::CellSearchMesh(
    Position low, Position hi, const std::array<uint32_t, 3>& shape,
    const std::vector<std::vector<std::shared_ptr<Cell>>>& elements)
    : shape_(shape),
      pitch_{0., 0., 0.},
      r_low(low),
      r_hi(hi),
      elements_(elements) {
  pitch_[0] = (r_hi.x() - r_low.x()) / static_cast<double>(shape_[0]);
  pitch_[1] = (r_hi.y() - r_low.y()) / static_cast<double>(shape_[1]);
  pitch_[2] = (r_hi.z() - r_low.z()) / static_cast<double>(shape_[2]);
}

std::array<int32_t, 3> CellSearchMesh::get_tile(const Position& r,
                                                const Direction& u) const {
  int32_t nx =
      static_cast<int32_t>(std::floor((r.x() - r_low.x()) / pitch_[0]));
  int32_t ny =
      static_cast<int32_t>(std::floor((r.y() - r_low.y()) / pitch_[1]));
  int32_t nz =
      static_cast<int32_t>(std::floor((r.z() - r_low.z()) / pitch_[2]));

  double x_c = r_low.x() + pitch_[0] * (static_cast<double>(nx) + 0.5);
  double y_c = r_low.y() + pitch_[1] * (static_cast<double>(ny) + 0.5);
  double z_c = r_low.z() + pitch_[2] * (static_cast<double>(nz) + 0.5);
  Position r_tile(x_c, y_c, z_c);

  // Bounds of tile. Must check all to see if we are on them, then use
  // direction to see if we are in or out of tile.
  double xl = r_tile.x() - pitch_[0] * 0.5;
  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT && u.x() < 0.) nx--;

  double xh = r_tile.x() + pitch_[0] * 0.5;
  if (std::abs(xh - r.x()) < SURFACE_COINCIDENT && u.x() >= 0.) nx++;

  double yl = r_tile.y() - pitch_[1] * 0.5;
  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT && u.y() < 0.) ny--;

  double yh = r_tile.y() + pitch_[1] * 0.5;
  if (std::abs(yh - r.y()) < SURFACE_COINCIDENT && u.y() >= 0.) ny++;

  double zl = r_tile.z() - pitch_[2] * 0.5;
  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT && u.z() < 0.) nz--;

  double zh = r_tile.z() + pitch_[2] * 0.5;
  if (std::abs(zh - r.z()) < SURFACE_COINCIDENT && u.z() >= 0.) nz++;

  return {nx, ny, nz};
}

bool CellSearchMesh::is_inside(const Position& r, const Direction& u) const {
  // Get index of each axis
  auto tile = get_tile(r, u);
  int nx = tile[0];
  int ny = tile[1];
  int nz = tile[2];

  if ((nx < 0 || nx >= static_cast<int>(shape_[0])) ||
      (ny < 0 || ny >= static_cast<int>(shape_[1])) ||
      (nz < 0 || nz >= static_cast<int>(shape_[2]))) {
    // Index is outside of mesh
    return false;
  }

  return true;
}

std::size_t CellSearchMesh::index(const Position& r, const Direction& u) const {
  auto tile = get_tile(r, u);

  if ((tile[0] < 0 || tile[0] >= static_cast<int>(shape_[0])) ||
      (tile[1] < 0 || tile[1] >= static_cast<int>(shape_[1])) ||
      (tile[2] < 0 || tile[2] >= static_cast<int>(shape_[2]))) {
    return -1;
  }

  std::size_t indx = shape_[2] * shape_[1] * static_cast<std::size_t>(tile[0]) +
                     shape_[1] * static_cast<std::size_t>(tile[1]) +
                     static_cast<std::size_t>(tile[2]);
  return indx;
}

const std::vector<std::shared_ptr<Cell>>& CellSearchMesh::index_cells(
    std::size_t i) const {
  return elements_[i];
}

Cell* CellSearchMesh::find_cell(const Position& r, const Direction& u,
                                int32_t on_surf) const {
  // First, we get our index
  std::size_t indx = this->index(r, u);

  if (indx >= this->size()) {
    // Apparently, we aren't in the mesh.
    return nullptr;
  }

  // We are in the mesh. Get the list of cells to check.
  const auto& cells = this->index_cells(indx);
  for (const auto& cell : cells) {
    if (cell->is_inside(r, u, on_surf)) return cell.get();
  }

  // Apparently we didn't find a cell we are in. This is odd. Probably
  // a bad mesh.
  return nullptr;
}

std::shared_ptr<CellSearchMesh> make_cell_search_mesh(const YAML::Node& node) {
  // Make sure the node is a map
  if (!node.IsMap()) {
    std::string mssg = "Cell search mesh entry must be a map.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Now get the shape
  std::array<uint32_t, 3> shape;
  if (!node["shape"] || !node["shape"].IsSequence() ||
      node["shape"].size() != 3) {
    std::string mssg =
        "No valid \"shape\" entry provided to cell search mesh entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  shape[0] = node["shape"][0].as<uint32_t>();
  shape[1] = node["shape"][1].as<uint32_t>();
  shape[2] = node["shape"][2].as<uint32_t>();

  // Now get the lower position
  if (!node["low"] || !node["low"].IsSequence() || node["low"].size() != 3) {
    std::string mssg =
        "No valid \"low\" entry provided to cell search mesh entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  double x_low = node["low"][0].as<double>();
  double y_low = node["low"][1].as<double>();
  double z_low = node["low"][2].as<double>();
  Position r_low(x_low, y_low, z_low);

  // Now get the upper position
  if (!node["hi"] || !node["hi"].IsSequence() || node["hi"].size() != 3) {
    std::string mssg =
        "No valid \"hi\" entry provided to cell search mesh entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  double x_hi = node["hi"][0].as<double>();
  double y_hi = node["hi"][1].as<double>();
  double z_hi = node["hi"][2].as<double>();
  Position r_hi(x_hi, y_hi, z_hi);

  // Now get all the elements
  if (!node["data"] || !node["data"].IsSequence() ||
      node["data"].size() != shape[0] * shape[1] * shape[2]) {
    std::string mssg =
        "No valid \"data\" entry provided to cell search mesh entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::vector<std::vector<uint32_t>> element_cell_ids =
      node["data"].as<std::vector<std::vector<uint32_t>>>();
  std::vector<std::vector<std::shared_ptr<Cell>>> elements(
      element_cell_ids.size(), std::vector<std::shared_ptr<Cell>>());

  // Go through all cells in all elements, and make sure they exist
  for (std::size_t ie = 0; ie < element_cell_ids.size(); ie++) {
    for (std::size_t ic = 0; ic < element_cell_ids[ie].size(); ic++) {
      if (cell_id_to_indx.find(element_cell_ids[ie][ic]) ==
          cell_id_to_indx.end()) {
        std::string mssg = "Cannot find cell with id " +
                           std::to_string(element_cell_ids[ie][ic]) +
                           " in cell-search-mesh element " +
                           std::to_string(ie) + ".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
      auto cell_indx = cell_id_to_indx[element_cell_ids[ie][ic]];
      elements[ie].push_back(geometry::cells[cell_indx]);
    }
  }

  return std::make_shared<CellSearchMesh>(r_low, r_hi, shape, elements);
}
