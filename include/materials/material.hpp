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
#ifndef MATERIAL_H
#define MATERIAL_H

#include <materials/nuclide.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

class Material;
extern std::map<uint32_t, std::shared_ptr<Material>> materials;

class Material : public std::enable_shared_from_this<Material> {
 public:
  enum class Fraction { Atom, Weight };
  enum class DensityUnits { g_cm3, a_bcm };

  struct NuclideInfo {
    std::string nuclide_key{};
    double awr{-1.};
    double atoms_bcm{-1.};
    double atom_fraction{-1.};
    double weight_fraction{-1.};

    NuclideInfo(const std::string& nuc_key) : nuclide_key(nuc_key) {}
  };

  struct Component {
    double atoms_bcm;  // atoms / (b * cm)
    std::shared_ptr<Nuclide> nuclide;
  };

  // This is used when making a ficticious material in noise vibrations
  Material() = default;

  // For initializing an empty material (i.e. plotting mode)
  Material(std::uint32_t id, const std::string& name = "");

  // For initializing a MG material.
  Material(std::shared_ptr<Nuclide> nuc, std::uint32_t id,
           const std::string& name = "");

  // For initializing a CE material.
  Material(double temperature, double density, DensityUnits dunits,
           Fraction frac, const std::vector<NuclideInfo>& nuclides,
           std::uint32_t id, const std::string& name = "");

  double max_energy() const {
    double max = 10000.;
    for (const auto& comp : components_) {
      if (comp.nuclide->max_energy() < max) {
        max = comp.nuclide->max_energy();
      }
    }

    return max;
  }

  double min_energy() const {
    double min = 0.;
    for (const auto& comp : components_) {
      if (comp.nuclide->min_energy() > min) {
        min = comp.nuclide->min_energy();
      }
    }

    return min;
  }

  // Used for building a custom material in noise vibrations
  const std::vector<Component>& components() const { return components_; }

  void add_component(const Component& comp) {
    if (comp.atoms_bcm <= 0.) {
      fatal_error("Cannot at component with negative atoms_bcm to a Material.");
    }

    components_.push_back(comp);
  }

  double temperature() const { return temperature_; }

  double density() const { return density_; }

  double atoms_per_bcm() const { return atoms_per_bcm_; }

  const std::string& name() const { return name_; }

  uint32_t id() const { return id_; }

  bool fissile() const {
    for (const auto& comp : components_) {
      if (comp.nuclide->fissile()) return true;
    }
    return false;
  }

 private:
  std::vector<NuclideInfo> nuclide_info_;
  std::vector<Component> components_;
  std::string name_;
  double temperature_;         // Temperature in Kelvin
  double density_{-1.};        // Density in g / cm^3
  double atoms_per_bcm_{-1.};  // atoms / (b * cm) sum of all components
  double avg_molar_mass_{-1.};
  std::uint32_t id_;

  double calc_avg_molar_mass(Fraction frac) const;

  friend void make_material(const YAML::Node& mat, bool plotting_mode);
  friend void fill_ce_material(const YAML::Node& mat,
                               std::shared_ptr<Material> material);
};

void make_material(const YAML::Node& mat, bool plotting_mode);

#endif
