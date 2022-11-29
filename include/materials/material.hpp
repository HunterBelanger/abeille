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

#include <map>
#include <materials/nuclide.hpp>
#include <memory>
#include <utils/rng.hpp>
#include <vector>

class Material;
extern std::map<uint32_t, std::shared_ptr<Material>> materials;

class Material : public std::enable_shared_from_this<Material> {
 public:
  Material() : composition_(), id_(), name_() {}

  struct Component {
    double concentration;
    std::shared_ptr<Nuclide> nuclide;
  };

  void add_component(double concentration, std::shared_ptr<Nuclide> nuclide) {
    composition_.push_back({concentration, nuclide});
  }

  double max_energy() const {
    double max = 10000.;
    for (const auto& comp : composition_) {
      if (comp.nuclide->max_energy() < max) {
        max = comp.nuclide->max_energy();
      }
    }

    return max;
  }

  double min_energy() const {
    double min = 0.;
    for (const auto& comp : composition_) {
      if (comp.nuclide->min_energy() > min) {
        min = comp.nuclide->min_energy();
      }
    }

    return min;
  }

  const std::vector<Component>& composition() const { return composition_; }
  const std::string& name() const { return name_; }
  uint32_t id() const { return id_; }
  bool fissile() const {
    for (const auto& comp : composition_) {
      if (comp.nuclide->fissile()) return true;
    }
    return false;
  }

 private:
  std::vector<Component> composition_;
  uint32_t id_;
  std::string name_;

  friend void make_material(const YAML::Node& mat, bool plotting_mode);
};

void make_material(const YAML::Node& mat, bool plotting_mode);

#endif
