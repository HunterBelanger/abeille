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
#include <geometry/cell_universe.hpp>
#include <geometry/geometry.hpp>
#include <utils/error.hpp>

CellUniverse::CellUniverse(std::vector<uint32_t> i_ind, uint32_t i_id,
                           std::string i_name)
    : Universe{i_id, i_name}, cell_indicies{i_ind} {}

Cell* CellUniverse::get_cell(Position r, Direction u, int32_t on_surf) const {
  // Go through each cell, and return the first one for which the
  // given position is inside the cell
  for (auto& indx : cell_indicies) {
    if (geometry::cells[indx]->is_inside(r, u, on_surf))
      return geometry::cells[indx].get();
  }
  // No cell found, particle is lost
  return nullptr;
}

Cell* CellUniverse::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                             Direction u, int32_t on_surf) const {
  // First push universe info onto the stack
  stack.push_back({GeoLilyPad::PadType::Universe, id_, r, {0, 0, 0}, false});

  // Go through each cell, and return the first one for which the
  // given position is inside the cell
  for (auto& indx : cell_indicies) {
    if (geometry::cells[indx]->is_inside(r, u, on_surf)) {
      auto cell_id = geometry::cells[indx]->id();

      // Save stack data for cell
      stack.push_back(
          {GeoLilyPad::PadType::Cell, cell_id, r, {0, 0, 0}, false});

      // Return cell
      return geometry::cells[indx].get();
    }
  }
  // No cell found, particle is lost
  return nullptr;
}

void make_cell_universe(YAML::Node uni_node) {
  // Get id
  uint32_t id;
  if (uni_node["id"] && uni_node["id"].IsScalar()) {
    id = uni_node["id"].as<uint32_t>();
  } else {
    std::string mssg = "Universe must have a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get name if present
  std::string name = "";
  if (uni_node["name"] && uni_node["name"].IsScalar()) {
    name = uni_node["name"].as<std::string>();
  }

  // Get cells
  std::vector<uint32_t> cells;
  if (uni_node["cells"] && uni_node["cells"].IsSequence()) {
    // Go through and check all cells
    for (size_t c = 0; c < uni_node["cells"].size(); c++) {
      uint32_t cell_id = uni_node["cells"][c].as<uint32_t>();

      // Get index for id
      uint32_t cell_indx = 0;
      if (cell_id_to_indx.find(cell_id) == cell_id_to_indx.end()) {
        std::string mssg = "Referenced cell id " + std::to_string(cell_id) +
                           " could not be found.";
        fatal_error(mssg, __FILE__, __LINE__);
      } else {
        cell_indx = static_cast<uint32_t>(cell_id_to_indx[cell_id]);
      }
      cells.push_back(cell_indx);
    }

  } else {
    std::string mssg = "Cell based universe must have a list of cells.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure universe id not already taken
  if (universe_id_to_indx.find(id) != universe_id_to_indx.end()) {
    std::string mssg =
        "Universe id " + std::to_string(id) + " appears multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make universe
  universe_id_to_indx[id] = geometry::universes.size();
  geometry::universes.push_back(
      std::make_shared<CellUniverse>(cells, id, name));
}
