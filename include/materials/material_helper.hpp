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
#ifndef MATERIAL_HELPER_H
#define MATERIAL_HELPER_H

#include <materials/material.hpp>
#include <utils/constants.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

class MaterialHelper {
 public:
  enum class BranchlessReaction { SCATTER, FISSION };

  MaterialHelper(Material* material, double E)
      : mat(nullptr), index(), xs(), E_(0.) {
    this->set_material(material, E);
  }

  void set_material(Material* material, double E) {
    if (mat != material) {
      mat = material;
      index.resize(mat->composition().size());
      xs.resize(mat->composition().size(), 0.);
      get_energy_index(E);  // This sets xs_set to false !
    } else if (E_ != E)
      get_energy_index(E);  // This sets xs_set to false !
  }

  double Et(double E, bool noise = false) {
    if (E_ != E) get_energy_index(E);

    double Et_ = 0;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      double comp_xs = composition(i).concentration *
                       composition(i).nuclide->total_xs(E_, index[i]);
      Et_ += comp_xs;
    }

    if (noise) {
      Et_ += Ew(E, noise);
    }

    return Et_;
  }

  double Ew(double E, bool noise = false) {
    if (E_ != E) get_energy_index(E);

    double Ew_ = 0.;

    if (noise) {
      double p_speed = speed(E, index[0]);
      Ew_ += settings::eta * settings::w_noise / p_speed;
    }

    return Ew_;
  }

  double Ea(double E) {
    if (E_ != E) get_energy_index(E);

    double Ea_ = 0;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      Ea_ += composition(i).concentration *
             (composition(i).nuclide->disappearance_xs(E_, index[i]) +
              composition(i).nuclide->fission_xs(E_, index[i]));
    }

    return Ea_;
  }

  double Es(double E) { return this->Et(E) - this->Ea(E); }

  double Ef(double E) {
    if (E_ != E) get_energy_index(E);

    double Ef_ = 0;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      if (composition(i).nuclide->fissile()) {
        Ef_ += composition(i).concentration *
               composition(i).nuclide->fission_xs(E_, index[i]);
      }
    }

    return Ef_;
  }

  double vEf(double E) {
    if (E_ != E) get_energy_index(E);

    double vEf_ = 0.;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      if (composition(i).nuclide->fissile()) {
        vEf_ += composition(i).concentration *
                composition(i).nuclide->nu_total(E_, index[i]) *
                composition(i).nuclide->fission_xs(E_, index[i]);
      }
    }

    return vEf_;
  }

  double Eelastic(double E) {
    if (E_ != E) get_energy_index(E);

    double Eelastic_ = 0;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      Eelastic_ += composition(i).concentration *
                   composition(i).nuclide->elastic_xs(E_, index[i]);
    }

    return Eelastic_;
  }

  double Emt(uint32_t mt, double E) {
    if (E_ != E) get_energy_index(E);

    double Emt_ = 0;

    for (size_t i = 0; i < mat->composition().size(); i++) {
      Emt_ += composition(i).concentration *
              composition(i).nuclide->reaction_xs(mt, E_, index[i]);
    }

    return Emt_;
  }

  std::pair<const Nuclide*, MicroXSs> sample_nuclide(double E, pcg32& rng,
                                                     bool noise = false) {
    if (E_ != E) {
      get_energy_index(E);
    }

    for (size_t i = 0; i < mat->composition().size(); i++) {
      double comp_xs = composition(i).concentration *
                       composition(i).nuclide->total_xs(E_, index[i]);
      xs[i] = comp_xs;
    }

    int comp_ind = RNG::discrete(rng, xs);
    size_t nuclide_index = index[comp_ind];
    const Nuclide* nuclide = composition(comp_ind).nuclide.get();

    // Get MicroXSs for use latter on
    MicroXSs xss;
    xss.energy = E_;
    xss.concentration = composition(comp_ind).concentration;
    xss.energy_index = nuclide_index;
    xss.total = nuclide->total_xs(xss.energy, xss.energy_index);
    xss.nu_total = nuclide->nu_total(xss.energy, xss.energy_index);
    xss.nu_delayed = nuclide->nu_delayed(xss.energy, xss.energy_index);
    xss.fission = nuclide->fission_xs(xss.energy, xss.energy_index);
    xss.absorption =
        nuclide->disappearance_xs(xss.energy, xss.energy_index) + xss.fission;

    if (noise) {
      xss.noise_copy =
          Ew(E, noise) /
          (xss.concentration * static_cast<double>(mat->composition().size()));
      xss.total += xss.noise_copy;
    }

    return {nuclide, xss};
  }

  std::pair<const Nuclide*, MicroXSs> sample_branchless_nuclide(
      double E, pcg32& rng, BranchlessReaction reaction) {
    if (E_ != E) {
      get_energy_index(E);
    }

    for (size_t i = 0; i < mat->composition().size(); i++) {
      double comp_xs = 0.;
      switch (reaction) {
        case BranchlessReaction::SCATTER:
          comp_xs = composition(i).concentration *
                    (composition(i).nuclide->total_xs(E_, index[i]) -
                     (composition(i).nuclide->disappearance_xs(E_, index[i]) +
                      composition(i).nuclide->fission_xs(E_, index[i])));
          break;
        case BranchlessReaction::FISSION:
          comp_xs = composition(i).concentration *
                    composition(i).nuclide->nu_total(E_, index[i]) *
                    composition(i).nuclide->fission_xs(E_, index[i]);
          break;
      }
      xs[i] = comp_xs;
    }

    int comp_ind = RNG::discrete(rng, xs);
    size_t nuclide_index = index[comp_ind];
    const Nuclide* nuclide = composition(comp_ind).nuclide.get();

    // Get MicroXSs for use latter on
    MicroXSs xss;
    xss.energy = E_;
    xss.concentration = composition(comp_ind).concentration;
    xss.energy_index = nuclide_index;
    xss.total = nuclide->total_xs(xss.energy, xss.energy_index);
    xss.nu_total = nuclide->nu_total(xss.energy, xss.energy_index);
    xss.nu_delayed = nuclide->nu_delayed(xss.energy, xss.energy_index);
    xss.fission = nuclide->fission_xs(xss.energy, xss.energy_index);
    xss.absorption =
        nuclide->disappearance_xs(xss.energy, xss.energy_index) + xss.fission;
    return {nuclide, xss};
  }

 private:
  Material* mat;
  std::vector<size_t> index;
  std::vector<double> xs;
  double E_;

  const Material::Component& composition(size_t i) {
    return mat->composition()[i];
  }

  void get_energy_index(double E) {
    E_ = E;

    // Go through all components and get energy index
    for (size_t i = 0; i < mat->composition().size(); i++) {
      index[i] = composition(i).nuclide->energy_grid_index(E_);
    }
  }

  // Function to get the speed of a particle in cm/s
  double speed(double E, std::size_t i) {
    // If we are in CE, it doesn't matter which nuclide we use
    // to get the speed. If we are in MG, then there should
    // only be 1 nuclide, so we always use component 0.
    return mat->composition()[0].nuclide->speed(E, i);
  }
};

#endif
