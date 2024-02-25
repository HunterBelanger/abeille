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
#include <materials/material.hpp>
#include <materials/mg_nuclide.hpp>
#include <materials/nuclide.hpp>
#include <plotting/slice_plot.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <memory>
#include <sstream>
#include <vector>

std::map<uint32_t, std::shared_ptr<Material>> materials;

Material::Material(std::uint32_t id, const std::string& name)
    : nuclide_info_(),
      components_(),
      name_(name),
      temperature_(293.6),
      density_(-1.),
      atoms_per_bcm_(1.),
      avg_molar_mass_(-1.),
      id_(id) {}

Material::Material(std::shared_ptr<Nuclide> nuc, std::uint32_t id,
                   const std::string& name)
    : nuclide_info_(),
      components_(),
      name_(name),
      temperature_(293.6),
      density_(-1.),
      atoms_per_bcm_(1.),
      avg_molar_mass_(-1.),
      id_(id) {
  // At the single nuclide component
  components_.push_back({1., nuc});
}

Material::Material(double temperature, double density, DensityUnits dunits,
                   Fraction frac, const std::vector<NuclideInfo>& nuclide_infos,
                   std::uint32_t id, const std::string& name)
    : nuclide_info_(nuclide_infos),
      components_(),
      name_(name),
      temperature_(temperature),
      density_(-1.),
      atoms_per_bcm_(1.),
      avg_molar_mass_(-1.),
      id_(id) {
  // First, check if temperatuer if positive
  if (temperature_ < 0.) {
    std::stringstream mssg;
    mssg << "Temperature provided to Material " << id_ << " was < 0.";
    fatal_error(mssg.str());
  }

  if (nuclide_info_.size() == 0 && density != 0.) {
    std::stringstream mssg;
    mssg << "Material " << id_ << " has no nuclides with a non-zero density.";
    fatal_error(mssg.str());
  } else if (nuclide_info_.size() != 0 && density == 0.) {
    std::stringstream mssg;
    mssg << "Material " << id_ << " has a density of zero with nuclides.";
    fatal_error(mssg.str());
  }

  // Now check the density
  if (density < 0.) {
    std::stringstream mssg;
    mssg << "Density provided to Material " << id_ << " was <= 0.";
    fatal_error(mssg.str());
  }

  // Check which density we were given
  if (dunits == DensityUnits::a_bcm) {
    atoms_per_bcm_ = density;
  } else {
    density_ = density;
  }

  if (nuclide_info_.size() == 0) {
    // Material is void
    density_ = 0.;
    avg_molar_mass_ = 0.;
    atoms_per_bcm_ = 0.;
  } else {
    // Material is not void

    // Make sure all awr and all fractions are positive
    for (const auto& nuc : nuclide_info_) {
      if (nuc.awr < 0.) {
        std::stringstream mssg;
        mssg << "Found negative awr for " << nuc.nuclide_key << " in material "
             << id_ << ".";
        fatal_error(mssg.str());
      }

      if ((frac == Fraction::Atom && nuc.atom_fraction < 0.) ||
          (frac == Fraction::Weight && nuc.weight_fraction < 0.)) {
        std::stringstream mssg;
        mssg << "Found negative fraction for " << nuc.nuclide_key
             << " in material " << id_ << ".";
        fatal_error(mssg.str());
      }
    }

    // Normalize given fractions
    double frac_sum = 0.;
    for (const auto& nuc : nuclide_info_) {
      if (frac == Fraction::Atom) {
        frac_sum += nuc.atom_fraction;
      } else {
        frac_sum += nuc.weight_fraction;
      }
    }
    for (auto& nuc : nuclide_info_) {
      if (frac == Fraction::Atom) {
        nuc.atom_fraction /= frac_sum;
      } else {
        nuc.weight_fraction /= frac_sum;
      }
    }

    // Now we need to get the average molar mass of the material before we can
    // get the other fraction.
    avg_molar_mass_ = this->calc_avg_molar_mass(frac);

    // Now we get the other fraction
    for (auto& nuc : nuclide_info_) {
      if (frac == Fraction::Atom) {
        nuc.weight_fraction =
            nuc.atom_fraction * nuc.awr * N_MASS_AMU / avg_molar_mass_;
      } else {
        nuc.atom_fraction =
            nuc.weight_fraction * avg_molar_mass_ / (nuc.awr * N_MASS_AMU);
      }
    }

    // Either denisty_ = -1. or atoms_per_bcm_ = -1. at this point. We need to
    // fix that one to the correct value
    if (density_ < 0. && atoms_per_bcm_ < 0.) {
      fatal_error("Something has gone terribly wrong.");
    } else if (density_ < 0.) {
      density_ = avg_molar_mass_ * atoms_per_bcm_ / N_AVAGADRO;
    } else {
      atoms_per_bcm_ = (density_ * N_AVAGADRO) / avg_molar_mass_;
    }

    // We can now set all of the atoms_bcm in each nuclide info !
    for (auto& nuc : nuclide_info_) {
      nuc.atoms_bcm = atoms_per_bcm_ * nuc.atom_fraction;
    }

    // At this point, all information in the class is set, and all of the
    // information in each NuclideInfo is set. We can now fill the components.
    for (const auto& nuc : nuclide_info_) {
      // Get CENuclidePacket from NDDirectory
      auto packet =
          settings::nd_directory->load_nuclide(nuc.nuclide_key, temperature_);

      if (!packet.nuclide_1 && !packet.nuclide_2) {
        // If we have no nuclides, this is bad
        // I don't think this can actually happen though. TODO
        std::stringstream mssg;
        mssg << "Nuclide " << nuc.nuclide_key << " in material ";
        mssg << id_ << " returned an empty CENuclidePacket.";
        fatal_error(mssg.str());
      }

      if (packet.nuclide_1) {
        const auto& temp_frac = packet.nuclide_1.value();

        // If this is a new nuclide, add it to the map of nuclides
        if (nuclides.find(temp_frac.nuclide->id()) == nuclides.end()) {
          nuclides[temp_frac.nuclide->id()] = temp_frac.nuclide;
        }

        // If this is a new zaid with a URR, add it to the zaids with urr set
        if (temp_frac.nuclide->has_urr() &&
            zaids_with_urr.find(temp_frac.nuclide->zaid()) == zaids_with_urr.end()) {
          zaids_with_urr.insert(temp_frac.nuclide->zaid());
        }

        components_.push_back(
            {nuc.atoms_bcm * temp_frac.fraction, temp_frac.nuclide});
      }

      if (packet.nuclide_2) {
        const auto& temp_frac = packet.nuclide_2.value();

        // If this is a new nuclide, add it to the map of nuclides
        if (nuclides.find(temp_frac.nuclide->id()) == nuclides.end()) {
          nuclides[temp_frac.nuclide->id()] = temp_frac.nuclide;
        }

        // If this is a new zaid with a URR, add it to the zaids with urr set
        if (temp_frac.nuclide->has_urr() &&
            zaids_with_urr.find(temp_frac.nuclide->zaid()) == zaids_with_urr.end()) {
          zaids_with_urr.insert(temp_frac.nuclide->zaid());
        }

        components_.push_back(
            {nuc.atoms_bcm * temp_frac.fraction, temp_frac.nuclide});
      }
    }
  }  // End if not void
}

double Material::calc_avg_molar_mass(Fraction frac) const {
  double avg_mm = 0.;

  for (const auto& nuc : nuclide_info_) {
    if (frac == Fraction::Atom) {
      avg_mm += nuc.atom_fraction * nuc.awr * N_MASS_AMU;
    } else {
      avg_mm += nuc.weight_fraction / (nuc.awr * N_MASS_AMU);
    }
  }

  if (frac == Fraction::Weight) {
    avg_mm = 1. / avg_mm;
  }

  return avg_mm;
}

void fill_mg_material(const YAML::Node& mat,
                      std::shared_ptr<Material> material) {
  std::string read_mssg = " Reading material " + material->name() + ".\n";
  Output::instance().write(read_mssg);

  // Can add nuclide to material
  std::shared_ptr<Nuclide> nuclide = make_mg_nuclide(mat, material->id());

  // Get the old name and id
  std::uint32_t id = material->id();
  std::string name = material->name();

  // Reconstruct material
  *material = Material(nuclide, id, name);
}

void fill_ce_material(const YAML::Node& mat,
                      std::shared_ptr<Material> material) {
  std::string read_mssg = " Reading material " + material->name() + ".\n";
  Output::instance().write(read_mssg);

  // First, get the temperature of the material
  double temp = 293.6;
  if (!mat["temperature"] || !mat["temperature"].IsScalar()) {
    std::stringstream mssg;
    mssg << "Material " << material->id() << " has no valid temperature entry.";
    fatal_error(mssg.str());
  }
  temp = mat["temperature"].as<double>();
  if (temp < 0.) {
    std::stringstream mssg;
    mssg << "Material " << material->id() << " has a negative temperature.";
    fatal_error(mssg.str());
  }

  // Now we get the fraction type. Assume atom fractions by default.
  Material::Fraction fractions = Material::Fraction::Atom;
  if (mat["fractions"] && !mat["fractions"].IsScalar()) {
    std::stringstream mssg;
    mssg << "Material " << material->id()
         << " has an ivalid \"fractions\" entry.";
    fatal_error(mssg.str());
  }
  if (mat["fractions"]) {
    std::string fractions_str = mat["fractions"].as<std::string>();

    if (fractions_str == "atoms") {
      fractions = Material::Fraction::Atom;
    } else if (fractions_str == "weight") {
      fractions = Material::Fraction::Weight;
    } else {
      std::stringstream mssg;
      mssg << "Material " << material->id()
           << " has invalid \"fractions\" entry.";
      fatal_error(mssg.str());
    }
  }

  // Now get density units str
  if (!mat["density-units"] || !mat["density-units"].IsScalar()) {
    std::stringstream mssg;
    mssg << "Material " << material->id()
         << " has invalid density-units entry.";
    fatal_error(mssg.str());
  }
  std::string density_units_str = mat["density-units"].as<std::string>();
  if (density_units_str != "g/cm3" && density_units_str != "atoms/b-cm" &&
      density_units_str != "sum") {
    std::stringstream mssg;
    mssg << "Material " << material->id() << " has invalid density units.";
    fatal_error(mssg.str());
  }

  Material::DensityUnits density_units = Material::DensityUnits::g_cm3;
  bool sum_for_density = false;
  if (density_units_str == "g/cm3") {
    density_units = Material::DensityUnits::g_cm3;
  } else if (density_units_str == "atoms/b-cm") {
    density_units = Material::DensityUnits::a_bcm;
  } else {
    sum_for_density = true;

    if (fractions == Material::Fraction::Atom) {
      density_units = Material::DensityUnits::a_bcm;
    } else {
      density_units = Material::DensityUnits::g_cm3;
    }
  }

  // Now get density if not sum
  double density = 0.;
  if (sum_for_density == false) {
    if (!mat["density"] || !mat["density"].IsScalar()) {
      std::stringstream mssg;
      mssg << "Material " << material->id() << " has invalid density entry.";
      fatal_error(mssg.str());
    }
    density = mat["density"].as<double>();

    if (density < 0.) {
      std::stringstream mssg;
      mssg << "Material " << material->id() << " has a nevative density.";
      fatal_error(mssg.str());
    }
  }

  // Finally, get the vectory of nuclides info
  if (!mat["composition"] || !mat["composition"].IsSequence()) {
    std::stringstream mssg;
    mssg << "Material " << material->id() << " has invalid composition list.";
    fatal_error(mssg.str());
  }
  auto comps = mat["composition"];
  std::vector<Material::NuclideInfo> nuclide_infos;

  for (std::size_t i = 0; i < comps.size(); i++) {
    // Make sure is a map !
    if (!comps[i].IsMap() || !comps[i]["nuclide"] ||
        !comps[i]["nuclide"].IsScalar() || !comps[i]["fraction"] ||
        !comps[i]["fraction"].IsScalar()) {
      std::stringstream mssg;
      mssg << "Entry index " << i << " in material ";
      mssg << material->id() << " is invalid.";
      fatal_error(mssg.str());
    }

    Material::NuclideInfo nuc_info(comps[i]["nuclide"].as<std::string>());
    nuc_info.awr = settings::nd_directory->awr(nuc_info.nuclide_key);
    if (fractions == Material::Fraction::Atom) {
      nuc_info.atom_fraction = comps[i]["fraction"].as<double>();
      if (sum_for_density) density += nuc_info.atom_fraction;

      if (nuc_info.atom_fraction < 0.) {
        std::stringstream mssg;
        mssg << "Entry " << i << " of material " << material->id() << " has a";
        mssg << " negative fraction.";
        fatal_error(mssg.str());
      }
    } else {
      nuc_info.weight_fraction = comps[i]["fraction"].as<double>();
      if (sum_for_density) density += nuc_info.weight_fraction;

      if (nuc_info.weight_fraction < 0.) {
        std::stringstream mssg;
        mssg << "Entry " << i << " of material " << material->id() << " has a";
        mssg << " negative fraction.";
        fatal_error(mssg.str());
      }
    }

    nuclide_infos.push_back(nuc_info);
  }

  // Get the name and id again
  std::uint32_t id = material->id();
  std::string name = material->name();

  // Reconstruct material
  *material = Material(temp, density, density_units, fractions, nuclide_infos,
                       id, name);
}

void make_material(const YAML::Node& mat, bool plotting_mode) {
  std::shared_ptr<Material> material = std::make_shared<Material>();

  // Get material id
  if (mat["id"] && mat["id"].IsScalar()) {
    material->id_ = mat["id"].as<uint32_t>();
  } else {
    fatal_error("Material is missing a valid id.");
  }

  // make sure material id isn't taken
  if (materials.find(material->id()) != materials.end()) {
    std::string mssg = "A material with id " + std::to_string(material->id()) +
                       " already exists.";
    fatal_error(mssg);
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
