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
#include <geometry/geometry.hpp>
#include <geometry/lattice.hpp>
#include <geometry/boundary.hpp>
#include <geometry/hex_lattice.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <utils/error.hpp>
#include <cstdint>

Lattice::Lattice(uint32_t i_id, std::string i_name)
    :  Universe{i_id, i_name},lattice_universes{}, outer_universe_index{-1}
    {
      this->has_boundary_conditions_ = false;

  if (this->outer_universe()) {
    this->has_boundary_conditions_ = this->outer_universe()->has_boundary_conditions();
  }
}

void Lattice::set_outisde_universe(int32_t univ) {
  outer_universe_index = univ;
}

Universe* Lattice::outer_universe() const {
  if (outer_universe_index < 0) return nullptr;

  return geometry::universes[static_cast<std::size_t>(outer_universe_index)]
      .get();
}

const Universe* Lattice::get_universe(std::size_t ind) const {
  if (lattice_universes[ind] < 0) return nullptr;

  std::size_t uni_indx = static_cast<std::size_t>(lattice_universes[ind]);

  return geometry::universes[uni_indx].get();
}



Boundary Lattice::get_boundary_condition(const Position& r,
                                                 const Direction& u,
                                                 int32_t on_surf) const {
  if (this->has_boundary_conditions()) {
    // Get lattice

    // If there is a boundary condition here, we know there is an outer
    // universe, so we are safe to dereference.
    return this->outer_universe()->get_boundary_condition(r, u, on_surf);
  }

  // No outer universe or no boundary conditions, so this goes out to infinity
  Boundary ret_bound(INF, -1, BoundaryType::Vacuum);
  ret_bound.token = 0;
  return ret_bound;
}

Boundary Lattice::lost_get_boundary(const Position& r,
                                            const Direction& u,
                                            int32_t on_surf) const {
  // Get lattice

  const bool is_inside = this->is_inside(r, u);
  if (this->outer_universe() && is_inside == false) {
    return this->outer_universe()->lost_get_boundary(r, u, on_surf);
  }

  std::array<int32_t, 3> tile = this->get_tile(r, u);
  Boundary ret_bound(this->distance_to_tile_boundary(r, u, tile), -1,
                     BoundaryType::Normal);
  ret_bound.token = 0;

  return ret_bound;
}

bool Lattice::contains_universe(uint32_t id) const {
  // Get lattice

  // First check outer
  if (this->outer_universe()) {
    if (this->outer_universe()->id() == id ||
        this->outer_universe()->contains_universe(id))
      return true;
  }

  // Now go through all other universes in lattice
  for (std::size_t ind = 0; ind < this->size(); ind++) {
    const Universe* uni = this->get_universe(ind);

    if (uni) {
      if (uni->id() == id || uni->contains_universe(id)) return true;
    }
  }

  return false;
}
/*
void LatticeUniverse::get_all_contained_cells()  
{
  Lattice* lat = geometry::lattices[lattice_index].get();
  for (std::size_t ind = 0; ind < lat->size(); ind++) {
     Universe* uni = lat->get_universe(ind);
    //if the cell has a universe in it
    if(uni)
      uni->get_all_contained_cells();
  }
}
*/

uint32_t Lattice::get_num_cell_instances(uint32_t cell_id) const {
   uint32_t instances = 0;

  // go through all tiles/universes
  for (std::size_t ind = 0; ind < this->size(); ind++) {
     const Universe* uni = this->get_universe(ind);    
     if(uni) {
      instances += uni->get_num_cell_instances(cell_id);
     }
  }

  if (this->outer_universe()) {
    instances += this->outer_universe()->get_num_cell_instances(cell_id);
  }

  return instances; 
}

std::set<uint32_t> Lattice::get_all_mat_cells() const {
  std::set<uint32_t> mat_cells;

  // go through all tiles/universes
  for (std::size_t ind = 0; ind < this->size(); ind++) {
    const Universe* uni = this->get_universe(ind);    
    if (uni) {
      auto uni_mat_cells = uni->get_all_mat_cells(); 
      mat_cells.insert(uni_mat_cells.begin(), uni_mat_cells.end());
    }
  }

  if (this->outer_universe()) {
      auto outer_mat_cells = this->outer_universe()->get_all_mat_cells(); 
      mat_cells.insert(outer_mat_cells.begin(), outer_mat_cells.end());
  }

  return mat_cells;
}

std::vector<std::map<const uint32_t, uint32_t>> Lattice::get_offset_map() const {
  
  std::vector<std::map<const uint32_t, uint32_t>> tile_id_offsets;
  tile_id_offsets.resize(this->size() + (this->outer_universe() ? 1 : 0));

  // Set of all material cells contained in this universe
  std::set<uint32_t> mat_cell_ids = this->get_all_mat_cells();

  // Set all offsets to be zero
  for (auto& map : tile_id_offsets) {
    for (uint32_t mat_cell_id : mat_cell_ids) {
      map[mat_cell_id] = 0;
    }
  }

  // Now go through and build all offsets
  for (std::size_t i = 1; i < this->size(); i++) {
    const Universe* uni = this->get_universe(i-1);

    for(uint32_t mat_cell_id : mat_cell_ids) {
      tile_id_offsets[i][mat_cell_id] = tile_id_offsets[i-1][mat_cell_id];

      if (uni) {
        tile_id_offsets[i][mat_cell_id] += uni->get_num_cell_instances(mat_cell_id);
      }
    }
  }

  // Also get offsets for outer universe
  if (this->outer_universe()) {
    auto uni = this->outer_universe();
    const std::size_t i = tile_id_offsets.size() - 1;

    for(uint32_t mat_cell_id : mat_cell_ids) {
      tile_id_offsets[i][mat_cell_id] = tile_id_offsets[i-1][mat_cell_id] + uni->get_num_cell_instances(mat_cell_id);
    }
  }

  return tile_id_offsets;
}

void make_lattice_universe(const YAML::Node& uni_node,
                           const YAML::Node& input) {
  // Get id
  uint32_t id = 0;
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

  // Get lattice
  uint32_t lat_id = 0;
  if (uni_node["lattice"] && uni_node["lattice"].IsScalar()) {
    lat_id = uni_node["lattice"].as<uint32_t>();
  } else {
    fatal_error("Lattice universe must have a valid lattice id.");
  }

  // See if lattice exists yet or not
  if (lattice_id_to_indx.find(lat_id) == lattice_id_to_indx.end()) {
    bool lattice_found = false;

    // Find lattice in input
    if (input["lattices"] && input["lattices"].IsSequence()) {
      // Iterate through lattices
      for (size_t l = 0; l < input["lattices"].size(); l++) {
        if (input["lattices"][l]["id"] &&
            input["lattices"][l]["id"].IsScalar()) {
          if (input["lattices"][l]["id"].as<uint32_t>() == lat_id) {
            // Lattice is a match, get type and then read it
            std::string type;
            if (input["lattices"][l]["type"] &&
                input["lattices"][l]["type"].IsScalar()) {
              type = input["lattices"][l]["type"].as<std::string>();
              if (type == "rectlinear") {
                make_rect_lattice(input["lattices"][l], input);
              } else if (type == "hexagonal") {
                make_hex_lattice(input["lattices"][l], input);
              } else {
                fatal_error("Unknown lattice type " + type + ".");
              }
            } else {
              fatal_error("Lattice instances must have a valid type.");
            }
            lattice_found = true;
            break;
          }
        } else {
          fatal_error("Lattice instances must have a valid id.");
        }
      }
      // If still not found, error
      if (lattice_found == false) {
        std::stringstream mssg;
        mssg << "Lattice with id " << lat_id << " could not be found.";
        fatal_error(mssg.str());
      }
    } else {
      fatal_error("Must have lattices in order to use a lattice universe.");
    }
  }

  // Make sure id isn't take
  if (lattice_id_to_indx.find(id) != lattice_id_to_indx.end()) {
    std::stringstream mssg;
    mssg << "Lattice id " << id << " appears multiple times.";
    fatal_error(mssg.str());
  }

  // Make lattice
  lattice_id_to_indx[id] = geometry::lattices.size();
  geometry::lattices.push_back(
      std::make_shared<Lattice>(id, name));
}
