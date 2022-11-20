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
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice_universe.hpp>
#include <geometry/rect_lattice.hpp>
#include <utils/error.hpp>

LatticeUniverse::LatticeUniverse(uint32_t i_lat, uint32_t i_id,
                                 std::string i_name)
    : Universe{i_id, i_name}, lattice_index{i_lat} {}

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

void make_lattice_universe(YAML::Node uni_node, YAML::Node input) {
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

  // Get lattice
  uint32_t lat_id;
  if (uni_node["lattice"] && uni_node["lattice"].IsScalar()) {
    lat_id = uni_node["lattice"].as<uint32_t>();
  } else {
    std::string mssg = "Lattice universe must have a valid lattice id.";
    fatal_error(mssg, __FILE__, __LINE__);
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
                std::string mssg = "Unknown lattice type " + type + ".";
                fatal_error(mssg, __FILE__, __LINE__);
              }
            } else {
              std::string mssg = "Lattice instances must have a valid type.";
              fatal_error(mssg, __FILE__, __LINE__);
            }
            lattice_found = true;
            break;
          }
        } else {
          std::string mssg = "Lattice instances must have a valid id.";
          fatal_error(mssg, __FILE__, __LINE__);
        }
      }
      // If still not found, error
      if (lattice_found == false) {
        std::string mssg = "Lattice with id " + std::to_string(lat_id) +
                           " could not be found.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    } else {
      std::string mssg =
          "Must have lattices in order to use a lattice universe.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Get lattice index
  uint32_t lat_indx = lattice_id_to_indx[lat_id];

  // Make sure id isn't take
  if (universe_id_to_indx.find(id) != universe_id_to_indx.end()) {
    std::string mssg =
        "Universe id " + std::to_string(id) + " appears multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make universe
  universe_id_to_indx[id] = geometry::universes.size();
  geometry::universes.push_back(
      std::make_shared<LatticeUniverse>(lat_indx, id, name));
}
