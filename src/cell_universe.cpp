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
#include <geometry/boundary.hpp>
#include <geometry/cell_universe.hpp>
#include <geometry/geometry.hpp>
#include <utils/error.hpp>

CellUniverse::CellUniverse(std::vector<uint32_t> i_ind, uint32_t i_id,
                           std::string i_name)
    : Universe{i_id, i_name}, cell_indicies{i_ind} {
  this->has_boundary_conditions_ = false;
  for (auto& indx : cell_indicies) {
    Cell* cell = geometry::cells[indx].get();
    if (cell->vacuum_or_reflective()) {
      this->has_boundary_conditions_ = true;
      break;
    }
  }
}

UniqueCell CellUniverse::get_cell(Position r, Direction u,
                                  int32_t on_surf) const {
  UniqueCell ucell;

  // Go through each cell, and return the first one for which the
  // given position is inside the cell
  for (std::size_t i = 0; i < cell_indicies.size(); i++) {
    const auto& indx = cell_indicies[i];

    if (geometry::cells[indx]->is_inside(r, u, on_surf)) {
      Cell* cell = geometry::cells[indx].get();

      if (cell->fill() == Cell::Fill::Material) {
        ucell.cell = cell;
        // This is a deep as it goes, so we set the ID here
        ucell.id = ucell.cell->id();
        ucell.instance += cell_offset_map[i].at(ucell.id);
        return ucell;
      }

      ucell = cell->universe()->get_cell(r, u, on_surf);
      ucell.instance += cell_offset_map[i].at(ucell.id);
      return ucell;
    }
  }
  // No cell found, particle is lost
  return ucell;
}

UniqueCell CellUniverse::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                                  Direction u, int32_t on_surf) const {
  // First push universe info onto the stack
  stack.push_back({GeoLilyPad::PadType::Universe, id_, r, {0, 0, 0}, false});

  UniqueCell ucell;

  // Go through each cell, and return the first one for which the
  // given position is inside the cell
  for (std::size_t i = 0; i < cell_indicies.size(); i++) {
    const auto& indx = cell_indicies[i];

    if (geometry::cells[indx]->is_inside(r, u, on_surf)) {
      auto cell_id = geometry::cells[indx]->id();

      // Save stack data for cell
      stack.push_back(
          {GeoLilyPad::PadType::Cell, cell_id, r, {0, 0, 0}, false});

      Cell* cell = geometry::cells[indx].get();

      if (cell->fill() == Cell::Fill::Material) {
        ucell.cell = cell;
        // This is a deep as it goes, so we set the ID here
        ucell.id = ucell.cell->id();
        ucell.instance += cell_offset_map[i].at(ucell.id);
        return ucell;
      }

      ucell = cell->universe()->get_cell(stack, r, u, on_surf);
      ucell.instance += cell_offset_map[i].at(ucell.id);
      return ucell;
    }
  }

  // No cell found, particle is lost
  return ucell;
}

Boundary CellUniverse::get_boundary_condition(const Position& r,
                                              const Direction& u,
                                              int32_t on_surf) const {
  double dist = INF;
  BoundaryType btype = BoundaryType::Vacuum;
  int surface_index = -1;
  int32_t token = 0;

  // Go through each cell, and check for boundary condition
  if (this->has_boundary_conditions()) {
    for (auto& indx : cell_indicies) {
      Cell* cell = geometry::cells[indx].get();

      // Only look at vacuum or reflective boundary conditions
      if (cell->vacuum_or_reflective() == false) continue;

      auto d_t = cell->distance_to_boundary_condition(r, u, on_surf);
      if (d_t.first < dist && std::abs(d_t.first - dist) > BOUNDRY_TOL) {
        double tmp_dist = d_t.first;
        int32_t tmp_token = std::abs(d_t.second);

        if (tmp_token) {
          token = tmp_token;
          dist = tmp_dist;
          surface_index = token - 1;
        } else {
          // Not an actuall surface
          continue;
        }

        btype = geometry::surfaces[static_cast<std::size_t>(surface_index)]
                    ->boundary();

        if (geometry::surfaces[static_cast<std::size_t>(surface_index)]->sign(
                r, u) < 0)
          token *= -1;
      }
    }
  }

  Boundary ret_bound(dist, surface_index, btype);
  ret_bound.token = token;
  return ret_bound;
}

bool CellUniverse::contains_universe(uint32_t id) const {
  for (auto& indx : cell_indicies) {
    Cell* cell = geometry::cells[indx].get();

    if (cell->fill() == Cell::Fill::Universe) {
      if (cell->universe()->id() == id ||
          cell->universe()->contains_universe(id))
        return true;
    }
  }

  return false;
}

Boundary CellUniverse::lost_get_boundary(const Position& r, const Direction& u,
                                         int32_t on_surf) const {
  double dist = INF;
  BoundaryType btype = BoundaryType::Vacuum;
  int surface_index = -1;
  int32_t token = 0;

  for (auto& indx : cell_indicies) {
    Cell* cell = geometry::cells[indx].get();

    // First check the boundary of the cell itself
    auto d_t = cell->distance_to_boundary(r, u, on_surf);
    if (d_t.first < dist && std::abs(d_t.first - dist) > BOUNDRY_TOL) {
      double tmp_dist = d_t.first;
      int32_t tmp_token = std::abs(d_t.second);

      if (tmp_token) {
        token = tmp_token;
        dist = tmp_dist;
        surface_index = token - 1;
      } else {
        // Not an actual surface
        continue;
      }

      btype = geometry::surfaces[static_cast<std::size_t>(surface_index)]
                  ->boundary();

      if (geometry::surfaces[static_cast<std::size_t>(surface_index)]->sign(
              r, u) < 0)
        token *= -1;
    }

    // Now check the universe in the cell, if there is one
    if (cell->fill() == Cell::Fill::Universe) {
      auto cell_uni_bound = cell->universe()->lost_get_boundary(r, u, on_surf);
      if (cell_uni_bound.distance < dist &&
          std::abs(cell_uni_bound.distance - dist) > BOUNDRY_TOL) {
        dist = cell_uni_bound.distance;
        surface_index = cell_uni_bound.surface_index;
        token = cell_uni_bound.token;
        btype = cell_uni_bound.boundary_type;
      }
    }
  }

  Boundary ret_bound(dist, surface_index, btype);
  ret_bound.token = token;
  return ret_bound;
}

std::set<uint32_t> CellUniverse::get_all_mat_cells() const {
  std::set<uint32_t> mat_cells;

  for (auto& indx : cell_indicies) {
    Cell* cell = geometry::cells[indx].get();

    if (cell->fill() == Cell::Fill::Material) {
      // Add the material cell to the set
      mat_cells.insert(cell->id());
    } else {
      // Get all material cells from the universe
      auto nested_mat_cells = cell->universe()->get_all_mat_cells();
      mat_cells.insert(nested_mat_cells.begin(), nested_mat_cells.end());
    }
  }

  return mat_cells;
}

uint32_t CellUniverse::get_num_cell_instances(uint32_t cell_id) const {
  uint32_t instances = 0;

  // go through all cells
  for (auto& indx : cell_indicies) {
    Cell* cell = geometry::cells[indx].get();

    if (cell->id() == cell_id) {
      instances++;
    }

    if (cell->fill() == Cell::Fill::Universe) {
      // check further universes
      instances += cell->universe()->get_num_cell_instances(cell_id);
    }
  }

  return instances;
}

void CellUniverse::make_offset_map() {
  cell_offset_map.resize(cell_indicies.size());
  // Set of all material cells contained in this universe
  std::set<uint32_t> mat_cell_ids = this->get_all_mat_cells();

  // Set offsets to all be zero initially
  for (auto& map : cell_offset_map) {
    for (uint32_t mat_cell_id : mat_cell_ids) {
      map[mat_cell_id] = 0;
    }
  }

  // Now go through and build all offsets
  for (std::size_t i = 1; i < cell_indicies.size(); i++) {
    for (uint32_t mat_cell_id : mat_cell_ids) {
      cell_offset_map[i][mat_cell_id] = cell_offset_map[i - 1][mat_cell_id];

      Cell* cell = geometry::cells[cell_indicies[i - 1]].get();
      if (cell->fill() == Cell::Fill::Universe) {
        cell_offset_map[i][mat_cell_id] +=
            cell->universe()->get_num_cell_instances(mat_cell_id);
      } else if (cell->id() == mat_cell_id) {
        // Currently, this isn't strictly necessary, as a material cell should
        // only appear in a universe once, so a non-zero offset would never be
        // used. This might become important, however, if cell transformations
        // are added one day.
        cell_offset_map[i][mat_cell_id] += 1;
      }
    }
  }
}

void make_cell_universe(const YAML::Node& uni_node) {
  // Get id
  uint32_t id;
  if (uni_node["id"] && uni_node["id"].IsScalar()) {
    id = uni_node["id"].as<uint32_t>();
  } else {
    fatal_error("Universe must have a valid id.");
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
        std::stringstream mssg;
        mssg << "Referenced cell id " << cell_id << " could not be found.";
        fatal_error(mssg.str());
      } else {
        cell_indx = static_cast<uint32_t>(cell_id_to_indx[cell_id]);
      }
      cells.push_back(cell_indx);
    }

  } else {
    fatal_error("Cell based universe must have a list of cells.");
  }

  // Make sure universe id not already taken
  if (universe_id_to_indx.find(id) != universe_id_to_indx.end()) {
    std::stringstream mssg;
    mssg << "Universe id " << id << " appears multiple times.";
    fatal_error(mssg.str());
  }

  // Make universe
  universe_id_to_indx[id] = geometry::universes.size();
  geometry::universes.push_back(
      std::make_shared<CellUniverse>(cells, id, name));
}
