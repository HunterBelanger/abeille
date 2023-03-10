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
#include <utils/rng.hpp>
#include <utils/settings.hpp>
#include <materials/nuclide.hpp>

#include <memory>
#include <sstream>
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

void fill_ce_material(const YAML::Node& mat,
                      std::shared_ptr<Material> material) {
  std::string read_mssg = " Reading material " + material->name() + ".\n";
  Output::instance()->write(read_mssg);

  // First, get the temperature of the material
  double temp = 293.6;
  if (!mat["temperature"] || !mat["temperature"].IsScalar()) {
    std::stringstream mssg;
    mssg << "Material " << material->name()
         << " has no valid temperature entry.";
    fatal_error(mssg.str(), __FILE__, __LINE__);
  }
  temp = mat["temperature"].as<double>();
  if (temp < 0.) {
    std::stringstream mssg;
    mssg << "Material " << material->name() << " has a negative temperature.";
    fatal_error(mssg.str(), __FILE__, __LINE__);
  }

  // Get YAML list of all components
  if (!mat["composition"] || !mat["composition"].IsSequence()) {
    std::stringstream mssg;
    mssg << "Material " << material->name() << " has invalid composition list.";
    fatal_error(mssg.str(), __FILE__, __LINE__);
  }
  auto comps = mat["composition"];

  for (std::size_t i = 0; i < comps.size(); i++) {
    // Make sure is a map !
    if (!comps[i].IsMap() || !comps[i]["nuclide"] ||
        !comps[i]["nuclide"].IsScalar() || !comps[i]["concentration"] ||
        !comps[i]["concentration"].IsScalar()) {
      std::stringstream mssg;
      mssg << "Entry index " << i << " in material ";
      mssg << material->name() << " is invalid.";
      fatal_error(mssg.str(), __FILE__, __LINE__);
    }

    const std::string nuclide_key = comps[i]["nuclide"].as<std::string>();
    const double concentration = comps[i]["concentration"].as<double>();

    // Make sure concentration is positive
    if (concentration < 0.) {
      std::stringstream mssg;
      mssg << "Entry index " << i << " in material ";
      mssg << material->name() << " has a negative concentration.";
      fatal_error(mssg.str(), __FILE__, __LINE__);
    }

    // Get CENuclidePacket from NDDirectory
    auto packet = settings::nd_directory->load_nuclide(nuclide_key, temp);

    // Add components. I don't think this can actually happen though. TODO
    if (!packet.nuclide_1 && !packet.nuclide_2) {
      // If we have no nuclides, this is bad
      std::stringstream mssg;
      mssg << "Entry index " << i << " in material ";
      mssg << material->name() << " returned an empty CENuclidePacket.";
      fatal_error(mssg.str(), __FILE__, __LINE__);
    }

    if (packet.nuclide_1) {
      const auto& nuc_frac = packet.nuclide_1.value();
      
      // If this is a new nuclide, add it to the map of nuclides
      if (nuclides.find(nuc_frac.nuclide->id()) == nuclides.end()) {
        nuclides[nuc_frac.nuclide->id()] = nuc_frac.nuclide;
      }

      // If this is a new zaid with a URR, add it to the zaids with urr set
      if (nuc_frac.nuclide->has_urr() &&
          zaids_with_urr.contains(nuc_frac.nuclide->zaid()) == false) {
        zaids_with_urr.insert(nuc_frac.nuclide->zaid());
      }

      material->add_component(nuc_frac.fraction * concentration,
                              nuc_frac.nuclide);
    }

    if (packet.nuclide_2) {
      const auto& nuc_frac = packet.nuclide_2.value();
      
      // If this is a new nuclide, add it to the map of nuclides
      if (nuclides.find(nuc_frac.nuclide->id()) == nuclides.end()) {
        nuclides[nuc_frac.nuclide->id()] = nuc_frac.nuclide;
      }

      // If this is a new zaid with a URR, add it to the zaids with urr set
      if (nuc_frac.nuclide->has_urr() &&
          zaids_with_urr.contains(nuc_frac.nuclide->zaid()) == false) {
        zaids_with_urr.insert(nuc_frac.nuclide->zaid());
      }

      material->add_component(nuc_frac.fraction * concentration,
                              nuc_frac.nuclide);
    }
  }
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
    switch (settings::energy_mode) {
      case settings::EnergyMode::MG:
        fill_mg_material(mat, material);
        break;

      case settings::EnergyMode::CE:
        fill_ce_material(mat, material);
        break;
    }
  }

  // Save material
  materials[material->id()] = material;
}
