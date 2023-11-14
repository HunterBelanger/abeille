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
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice.hpp>
#include <geometry/lattice_universe.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <utils/error.hpp>

LatticeUniverse::LatticeUniverse(uint32_t i_lat, uint32_t i_id,
                                 std::string i_name)
    : Universe{i_id, i_name}, lattice_index{i_lat} {
  this->has_boundary_conditions_ = false;

  if (geometry::lattices[lattice_index]->outer_universe()) {
    this->has_boundary_conditions_ = geometry::lattices[lattice_index]
                                         ->outer_universe()
                                         ->has_boundary_conditions();
  }
}

Cell* LatticeUniverse::get_cell(Position r, Direction u,
                                int32_t on_surf) const {
  return geometry::lattices[lattice_index]->get_cell(r, u, on_surf);
}

Cell* LatticeUniverse::get_cell(std::vector<GeoLilyPad>& stack, Position r,
                                Direction u, int32_t on_surf) const {
  // First push universe info onto the stack
  stack.push_back({GeoLilyPad::PadType::Universe, id_, r, {0, 0, 0}, false});

  return geometry::lattices[lattice_index]->get_cell(stack, r, u, on_surf);
}

Boundary LatticeUniverse::get_boundary_condition(const Position& r,
                                                 const Direction& u,
                                                 int32_t on_surf) const {
  if (this->has_boundary_conditions()) {
    // Get lattice
    Lattice* lat = geometry::lattices[lattice_index].get();

    // If there is a boundary condition here, we know there is an outer
    // universe, so we are safe to dereference.
    return lat->outer_universe()->get_boundary_condition(r, u, on_surf);
  }

  // No outer universe or no boundary conditions, so this goes out to infinity
  Boundary ret_bound(INF, -1, BoundaryType::Vacuum);
  ret_bound.token = 0;
  return ret_bound;
}

Boundary LatticeUniverse::lost_get_boundary(const Position& r,
                                            const Direction& u,
                                            int32_t on_surf) const {
  // Get lattice
  Lattice* lat = geometry::lattices[lattice_index].get();

  const bool is_inside = lat->is_inside(r, u);
  if (lat->outer_universe() && is_inside == false) {
    return lat->outer_universe()->lost_get_boundary(r, u, on_surf);
  }

  std::array<int32_t, 3> tile = lat->get_tile(r, u);
  Boundary ret_bound(lat->distance_to_tile_boundary(r, u, tile), -1,
                     BoundaryType::Normal);
  ret_bound.token = 0;

  return ret_bound;
}

bool LatticeUniverse::contains_universe(uint32_t id) const {
  // Get lattice
  Lattice* lat = geometry::lattices[lattice_index].get();

  // First check outer
  if (lat->outer_universe()) {
    if (lat->outer_universe()->id() == id ||
        lat->outer_universe()->contains_universe(id))
      return true;
  }

  // Now go through all other universes in lattice
  for (std::size_t ind = 0; ind < lat->size(); ind++) {
    const Universe* uni = lat->get_universe(ind);

    if (uni) {
      if (uni->id() == id || uni->contains_universe(id)) return true;
    }
  }

  return false;
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

  // Get lattice index
  uint32_t lat_indx = static_cast<uint32_t>(lattice_id_to_indx[lat_id]);

  // Make sure id isn't take
  if (universe_id_to_indx.find(id) != universe_id_to_indx.end()) {
    std::stringstream mssg;
    mssg << "Universe id " << id << " appears multiple times.";
    fatal_error(mssg.str());
  }

  // Make universe
  universe_id_to_indx[id] = geometry::universes.size();
  geometry::universes.push_back(
      std::make_shared<LatticeUniverse>(lat_indx, id, name));
}
