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
#include <geometry/cell.hpp>
#include <geometry/geometry.hpp>
#include <geometry/surfaces/surface.hpp>
#include <plotting/plotter.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/parser.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <sstream>

Cell::Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Material> material,
           uint32_t i_id, std::string i_name)
    : fill_{Fill::Material},
      rpn{i_rpn},
      id_{i_id},
      name_{i_name},
      material_{material},
      material_raw_{material.get()},
      universe_{nullptr},
      universe_raw_{nullptr} {
  // Check if simple or not
  simplify();

  // Check for vacuum or reflective boundary conditions
  check_for_bc();
}

Cell::Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Universe> universe,
           uint32_t i_id, std::string i_name)
    : fill_{Fill::Universe},
      rpn{i_rpn},
      id_{i_id},
      name_{i_name},
      material_{nullptr},
      material_raw_{nullptr},
      universe_{universe},
      universe_raw_{universe.get()} {
  // Check if simple or not
  simplify();

  // Check for vacuum or reflective boundary conditions
  check_for_bc();
}

bool Cell::is_inside(const Position& r, const Direction& u,
                     int32_t on_surf) const {
  if (simple)
    return is_inside_simple(r, u, on_surf);
  else
    return is_inside_complex(r, u, on_surf);
}

std::pair<double, int32_t> Cell::distance_to_boundary(const Position& r,
                                                      const Direction& u,
                                                      int32_t on_surf) const {
  double min_dist = INF;
  int32_t i_surf{0};

  for (int32_t token : rpn) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP::UNIN) continue;

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident = std::abs(token) == std::abs(on_surf);
    double d =
        geometry::surfaces[static_cast<std::size_t>(abs(token) - 1)]->distance(
            r, u, coincident);

    // Check if this distance is the new minimum.
    if (d < min_dist) {
      if (std::abs(d - min_dist) / min_dist >= 1e-14) {
        min_dist = d;
        i_surf = -token;
      }
    }
  }

  return {min_dist, i_surf};
}

std::pair<double, int32_t> Cell::distance_to_boundary_condition(
    const Position& r, const Direction& u, int32_t on_surf) const {
  if (this->vacuum_or_reflective_ == false) return {INF, 0};

  double min_dist = INF;
  int32_t i_surf{0};

  for (int32_t token : rpn) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP::UNIN) continue;

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident = std::abs(token) == std::abs(on_surf);

    Surface* surf =
        geometry::surfaces[static_cast<std::size_t>(abs(token) - 1)].get();

    // Ignore surfaces which aren't vacuum or reflective
    if (surf->boundary() == BoundaryType::Normal) continue;
    ;

    double d = surf->distance(r, u, coincident);

    // Check if this distance is the new minimum.
    if (d < min_dist) {
      if (std::abs(d - min_dist) / min_dist >= 1e-14) {
        min_dist = d;
        i_surf = -token;
      }
    }
  }

  return {min_dist, i_surf};
}

bool Cell::is_inside_simple(const Position& r, const Direction& u,
                            int32_t on_surf) const {
  for (const int32_t& token : rpn) {
    if (token == on_surf) {
    } else if (-token == on_surf)
      return false;
    else {
      int sign =
          geometry::surfaces[static_cast<std::size_t>(std::abs(token) - 1)]
              ->sign(r, u);
      if ((sign > 0 && token < 0) || (sign < 0 && token > 0)) return false;
    }
  }
  return true;
}

bool Cell::is_inside_complex(const Position& r, const Direction& u,
                             int32_t on_surf) const {
  std::vector<bool> stck(rpn.size());
  int i_stck = -1;

  for (int32_t token : rpn) {
    if (token == OP::UNIN) {
      stck[static_cast<std::size_t>(i_stck - 1)] =
          stck[static_cast<std::size_t>(i_stck - 1)] ||
          stck[static_cast<std::size_t>(i_stck)];
      i_stck--;
    } else if (token == OP::INTR) {
      stck[static_cast<std::size_t>(i_stck - 1)] =
          stck[static_cast<std::size_t>(i_stck - 1)] &&
          stck[static_cast<std::size_t>(i_stck)];
      i_stck--;
    } else if (token == OP::COMP) {
      stck[static_cast<std::size_t>(i_stck)] =
          !stck[static_cast<std::size_t>(i_stck)];
    } else {
      i_stck++;
      if (token == on_surf) {
        stck[static_cast<std::size_t>(i_stck)] = true;
      } else if (-token == on_surf) {
        stck[static_cast<std::size_t>(i_stck)] = false;
      } else {
        int sign =
            geometry::surfaces[static_cast<std::size_t>(std::abs(token) - 1)]
                ->sign(r, u);
        if ((sign > 0 && token > 0) || (sign < 0 && token < 0)) {
          stck[static_cast<std::size_t>(i_stck)] = true;
        } else
          stck[static_cast<std::size_t>(i_stck)] = false;
      }
    }
  }

  if (i_stck == 0)
    return stck[static_cast<std::size_t>(i_stck)];
  else
    return true;
}

uint32_t Cell::id() const { return id_; }

const std::string& Cell::name() const { return name_; }

void Cell::check_for_bc() {
  // Check for vacuum or reflective boundary conditions
  for (int32_t token : rpn) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP::UNIN) continue;

    Surface* surf =
        geometry::surfaces[static_cast<std::size_t>(abs(token) - 1)].get();
    if (surf->boundary() == BoundaryType::Vacuum ||
        surf->boundary() == BoundaryType::Reflective) {
      vacuum_or_reflective_ = true;
    }
  }
}

void Cell::simplify() {
  // Check if simple or not
  simple = true;
  for (const auto& el : rpn) {
    if (el == OP::COMP || el == OP::UNIN) {
      simple = false;
      break;
    }
  }

  // If simple, remove un-needed operators
  if (simple) {
    size_t i0 = 0;
    size_t i1 = 0;
    while (i1 < rpn.size()) {
      if (rpn[i1] < OP::UNIN) {
        rpn[i0] = rpn[i1];
        i0++;
      }
      i1++;
    }
    rpn.resize(i0);
  }
  rpn.shrink_to_fit();
}

//============================================================================
// Non-Member functions
std::vector<int32_t> infix_to_rpn(const std::vector<int32_t>& infix) {
  std::vector<int32_t> rpn;
  std::vector<int32_t> stack;

  for (const auto& token : infix) {
    if (token < OP::UNIN) {
      rpn.push_back(token);
    } else if (token < OP::R_PAR) {
      while (stack.size() > 0) {
        int32_t op = stack.back();

        if (op < OP::R_PAR && ((token == OP::COMP && token < op) ||
                               (token != OP::COMP && token <= op))) {
          rpn.push_back(op);
          stack.pop_back();
        } else {
          break;
        }
      }

      stack.push_back(token);

    } else if (token == OP::L_PAR) {
      stack.push_back(token);
    } else {
      for (auto it = stack.rbegin(); *it != OP::L_PAR; it++) {
        if (it == stack.rend()) {
          fatal_error("Mismatched parentheses in cell region definition.");
        }
        rpn.push_back(stack.back());
        stack.pop_back();
      }
      stack.pop_back();
    }
  }

  while (stack.size() > 0) {
    int32_t op = stack.back();

    if (op >= OP::R_PAR) {
      // Thow error, parenth mis-match
      fatal_error("Mismatched parentheses in cell region definition.");
    }

    rpn.push_back(stack.back());
    stack.pop_back();
  }

  return rpn;
}

void make_cell(const YAML::Node& cell_node, const YAML::Node& input) {
  // Get region string
  std::string region_str;
  if (cell_node["region"]) {
    region_str = cell_node["region"].as<std::string>();
  } else {
    fatal_error("Cell missing region definition.");
  }

  // Parse region stirng
  std::vector<int32_t> region;
  std::string temp = "";
  for (size_t k = 0; k < region_str.size(); k++) {
    char c = region_str[k];
    if (c == '&' || c == '(' || c == ')' || c == 'U' || c == '~') {
      // Make sure temp not empty
      if (temp.size() > 0) {
        int32_t signed_id = std::stoi(temp);
        int32_t indx = static_cast<int32_t>(
            surface_id_to_indx[static_cast<uint32_t>(std::abs(signed_id))]);
        indx += 1;  // This is due to 1 off indexing of surfaces for use of the
                    // sign of tokens
        if (signed_id < 0) indx *= -1;
        region.push_back(indx);
        temp = "";
      }

      if (c == '&')
        region.push_back(OP::INTR);
      else if (c == '(')
        region.push_back(OP::L_PAR);
      else if (c == ')')
        region.push_back(OP::R_PAR);
      else if (c == 'U')
        region.push_back(OP::UNIN);
      else if (c == '~')
        region.push_back(OP::COMP);

    } else if ((c == '+') || (c == '-') || (c == '0') || (c == '1') ||
               (c == '2') || (c == '3') || (c == '4') || (c == '5') ||
               (c == '5') || (c == '6') || (c == '7') || (c == '8') ||
               (c == '9')) {
      temp += c;
    } else if (c != ' ') {
      // Invalid char
      fatal_error("Invalid character in cell region definition.");
    }
  }
  if (temp.size() > 0) {
    int32_t signed_id = std::stoi(temp);
    int32_t indx = static_cast<int32_t>(
        surface_id_to_indx[static_cast<uint32_t>(std::abs(signed_id))]);
    indx += 1;
    if (signed_id < 0) indx *= -1;
    region.push_back(indx);
    temp = "";
  }
  // Change from infix to rpn
  region = infix_to_rpn(region);

  // Get id
  uint32_t id = 0;
  if (cell_node["id"] && cell_node["id"].IsScalar()) {
    id = cell_node["id"].as<uint32_t>();
  } else {
    fatal_error("Cell is missing id.");
  }

  // Get name
  std::string name;
  if (cell_node["name"]) {
    name = cell_node["name"].as<std::string>();
  } else
    name = "";

  // Check if there is a material and cell definition
  if (cell_node["material"] && cell_node["universe"]) {
    std::stringstream mssg;
    mssg << "Cell with id " << id << " has both a material and a universe.";
    fatal_error(mssg.str());
  }

  if (!cell_node["material"] && !cell_node["universe"]) {
    std::stringstream mssg;
    mssg << "Cell with id " << id << " neither a material nor a universe.";
    fatal_error(mssg.str());
  }

  // Initial nullptr, will be filled with either universe or material cell
  std::shared_ptr<Cell> cell_pntr{nullptr};

  if (cell_node["material"]) {
    // Get material id
    uint32_t mat_id = 0;
    if (cell_node["material"] && cell_node["material"].IsScalar()) {
      mat_id = cell_node["material"].as<uint32_t>();
    } else {
      std::stringstream mssg;
      mssg << "Cell " << id << " has an invalid material definition.";
      fatal_error(mssg.str());
    }
    std::shared_ptr<Material> material = nullptr;
    // Make sure material exists
    if (materials.find(mat_id) == materials.end()) {
      std::stringstream mssg;
      mssg << "Could not find material with ID " << mat_id << ".";
      fatal_error(mssg.str());
    }
    material = materials[mat_id];

    cell_pntr = std::make_shared<Cell>(region, material, id, name);
  } else {
    // We must have a universe then
    // Get universe id
    uint32_t uni_id = 0;
    if (cell_node["universe"] && cell_node["universe"].IsScalar()) {
      uni_id = cell_node["universe"].as<uint32_t>();
    } else {
      std::stringstream mssg;
      mssg << "Cell " << id << " has an invalid material definition.";
      fatal_error(mssg.str());
    }

    // Make sure universe exists
    if (universe_id_to_indx.find(uni_id) == universe_id_to_indx.end()) {
      find_universe(input, uni_id);
    }

    auto uni_indx = universe_id_to_indx[uni_id];
    std::shared_ptr<Universe> universe = geometry::universes[uni_indx];

    cell_pntr = std::make_shared<Cell>(region, universe, id, name);
  }

  // Add cell ID to map of surface indicies
  if (cell_id_to_indx.find(cell_pntr->id()) != cell_id_to_indx.end()) {
    // ID already exists
    std::stringstream mssg;
    mssg << "The cell ID " << cell_pntr->id() << " appears multiple times.";
    fatal_error(mssg.str());
  } else {
    cell_id_to_indx[cell_pntr->id()] = geometry::cells.size();
    geometry::cells.push_back(cell_pntr);

    // Generate a random color for the cell
    uint8_t r = static_cast<uint8_t>(255.0 * settings::rng());
    uint8_t g = static_cast<uint8_t>(255.0 * settings::rng());
    uint8_t b = static_cast<uint8_t>(255.0 * settings::rng());
    plotter::Pixel cell_color(r, g, b);
    plotter::cell_id_to_color[cell_pntr->id()] = cell_color;
  }
}
