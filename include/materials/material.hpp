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
#ifndef MATERIAL_H
#define MATERIAL_H

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <materials/nuclide.hpp>
#include <memory>
#include <sstream>
#include <string>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
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
