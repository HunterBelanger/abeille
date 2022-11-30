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
#include <geometry/cell.hpp>
#include <plotting/plotter.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

Cell::Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Material> material,
           uint32_t i_id, std::string i_name)
    : rpn{},
      id_{i_id},
      name_{i_name},
      material_{material},
      material_raw_{material.get()} {
  // Check if simple or not
  simple = true;
  for (const auto& el : i_rpn) {
    if (el == OP::COMP || el == OP::UNIN) {
      simple = false;
      break;
    }
  }

  rpn = i_rpn;

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
    double d = geometry::surfaces[abs(token) - 1]->distance(r, u, coincident);

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
      int sign = geometry::surfaces[std::abs(token) - 1]->sign(r, u);
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
      stck[i_stck - 1] = stck[i_stck - 1] || stck[i_stck];
      i_stck--;
    } else if (token == OP::INTR) {
      stck[i_stck - 1] = stck[i_stck - 1] && stck[i_stck];
      i_stck--;
    } else if (token == OP::COMP) {
      stck[i_stck] = !stck[i_stck];
    } else {
      i_stck++;
      if (token == on_surf) {
        stck[i_stck] = true;
      } else if (-token == on_surf) {
        stck[i_stck] = false;
      } else {
        int sign = geometry::surfaces[std::abs(token) - 1]->sign(r, u);
        if ((sign > 0 && token > 0) || (sign < 0 && token < 0)) {
          stck[i_stck] = true;
        } else
          stck[i_stck] = false;
      }
    }
  }

  if (i_stck == 0)
    return stck[i_stck];
  else
    return true;
}

uint32_t Cell::id() const { return id_; }

const std::string& Cell::name() const { return name_; }

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
          std::string mssg =
              "Mismatched parentheses in cell region definition.";
          fatal_error(mssg, __FILE__, __LINE__);
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
      std::string mssg = "Mismatched parentheses in cell region definition.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    rpn.push_back(stack.back());
    stack.pop_back();
  }

  return rpn;
}

void make_cell(YAML::Node cell_node) {
  // Get region string
  std::string region_str;
  if (cell_node["region"]) {
    region_str = cell_node["region"].as<std::string>();
  } else {
    std::string mssg = "Cell missing region definition.";
    fatal_error(mssg, __FILE__, __LINE__);
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
        int32_t indx =
            static_cast<int32_t>(surface_id_to_indx[std::abs(signed_id)]);
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
      std::string mssg = "Invalid character in cell region definition.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }
  if (temp.size() > 0) {
    int32_t signed_id = std::stoi(temp);
    int32_t indx =
        static_cast<int32_t>(surface_id_to_indx[std::abs(signed_id)]);
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
    std::string mssg = "Cell is missing id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get material id
  uint32_t mat_id = 0;
  if (cell_node["material"] && cell_node["material"].IsScalar()) {
    mat_id = cell_node["material"].as<uint32_t>();
  } else {
    std::string mssg = "Cell is missing material id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::shared_ptr<Material> material = nullptr;
  // Make sure material exists
  if (materials.find(mat_id) == materials.end()) {
    std::string mssg =
        "Could not find material with ID " + std::to_string(mat_id) + ".";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  material = materials[mat_id];

  // Get name
  std::string name;
  if (cell_node["name"]) {
    name = cell_node["name"].as<std::string>();
  } else
    name = "";

  std::shared_ptr<Cell> cell_pntr =
      std::make_shared<Cell>(region, material, id, name);

  // Add cell ID to map of surface indicies
  if (cell_id_to_indx.find(cell_pntr->id()) != cell_id_to_indx.end()) {
    // ID already exists
    std::string mssg = "The cell ID " + std::to_string(cell_pntr->id()) +
                       " appears multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else {
    cell_id_to_indx[cell_pntr->id()] = geometry::cells.size();
    geometry::cells.push_back(cell_pntr);

    // Generate a random color for the cell
    uint8_t r = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    uint8_t g = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    uint8_t b = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    plotter::Pixel cell_color(r, g, b);
    plotter::cell_id_to_color[cell_pntr->id()] = cell_color;
  }
}
