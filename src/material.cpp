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
#include <materials/material.hpp>
#include <materials/mg_nuclide.hpp>
#include <plotting/slice_plot.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <utils/rng.hpp>
#include <memory>
#include <vector>

std::map<uint32_t, std::shared_ptr<Material>> materials;

void fill_mg_material(const YAML::Node& mat,
                      std::shared_ptr<Material> material) {
  std::string read_mssg = " Reading material " + material->name() + ".\n";
  Output::instance()->write(read_mssg);

  // Can add nuclide to material
  std::shared_ptr<Nuclide> nuclide = make_mg_nuclide(mat, material->id());

  material->add_component(1., nuclide);
}

void make_material(const YAML::Node& mat, bool plotting_mode) {
  std::shared_ptr<Material> material = std::make_shared<Material>();

  // Get material id
  if (mat["id"] && mat["id"].IsScalar()) {
    material->id_ = mat["id"].as<uint32_t>();
  } else {
    std::string mssg = "Material is missing a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // make sure material id isn't taken
  if (materials.find(material->id()) != materials.end()) {
    std::string mssg = "A material with id " + std::to_string(material->id()) +
                       " already exists.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get material name
  if (mat["name"] && mat["name"].IsScalar()) {
    material->name_ = mat["name"].as<std::string>();
  }

  // Get the color for the material if one is provided. Only needed for
  // plotting purposes
  plotter::Pixel mat_color;
  if (mat["color"] && mat["color"].IsSequence() && mat["color"].size() == 3) {
    uint8_t R = mat["color"][0].as<uint8_t>();
    uint8_t G = mat["color"][1].as<uint8_t>();
    uint8_t B = mat["color"][2].as<uint8_t>();
    mat_color = plotter::Pixel(R, G, B);
  } else {
    // Generate a random color for the material
    uint8_t r = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    uint8_t g = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    uint8_t b = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
    mat_color = plotter::Pixel(r, g, b);
  }
  plotter::material_id_to_color[material->id()] = mat_color;

  if (!plotting_mode) {
    if (settings::energy_mode == settings::EnergyMode::MG) {
      fill_mg_material(mat, material);
    } else if (settings::energy_mode == settings::EnergyMode::CE) {
      std::string mssg = "Continuous Energy Mode not supported.";
      fatal_error(mssg, __FILE__, __LINE__);
    } else {
      std::string mssg = "Unknown energy mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Save material
  materials[material->id()] = material;
}
