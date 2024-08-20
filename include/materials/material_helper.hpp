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
#ifndef MATERIAL_HELPER_H
#define MATERIAL_HELPER_H

#include <materials/material.hpp>
#include <materials/nuclide.hpp>
#include <utils/constants.hpp>
#include <utils/noise_parameters.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <boost/unordered/unordered_flat_map.hpp>

#include <cstdint>
#include <optional>

class MaterialHelper {
 public:
  enum class BranchlessReaction { SCATTER, FISSION };

  MaterialHelper(Material* material, double E,
                 std::optional<NoiseParameters> noise_params = std::nullopt);

  std::optional<NoiseParameters> noise_params() const { return noise_params_; }
  void set_noise_params(std::optional<NoiseParameters> np) {
    noise_params_ = np;
  }

  void set_material(Material* material, double E) {
    mat = material;
    this->set_energy(E);
  }

  void set_urr_rand_vals(RNG& rng) {
    for (auto& rand : zaid_to_urr_rand_) rand.second = rng();
    this->clear_xs();
  }

  void set_urr_rand_vals(
      const boost::unordered_flat_map<uint32_t, std::optional<double>>& vals) {
    zaid_to_urr_rand_ = vals;
    this->clear_xs();
  }

  const boost::unordered_flat_map<uint32_t, std::optional<double>>&
  urr_rand_vals() const {
    return zaid_to_urr_rand_;
  }

  void clear_urr_rand_vals() {
    for (auto& rand : zaid_to_urr_rand_) rand.second = std::nullopt;
    this->clear_xs();
  }

  double Et(double E) {
    this->set_energy(E);

    double Et_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Et_ += comp.atoms_bcm * micro_xs.total;
    }

    // Get noise copy component (0 if noise_params_ is empty)
    Et_ += this->Ew(E);

    return Et_;
  }

  double Ew(double E) {
    this->set_energy(E);

    double Ew_ = 0.;

    if (noise_params_ && mat->components().size() > 0) {
      std::size_t energy_index = 0;

      // Needing to check MG/CE here isn't exactly elegant, but I don't
      // have a better solution for now...
      if (settings::energy_mode == settings::EnergyMode::MG) {
        energy_index = mat->components()[0].nuclide->energy_grid_index(E);
      }

      double p_speed = speed(E, energy_index);
      Ew_ += noise_params_->eta * noise_params_->omega / p_speed;
    }

    return Ew_;
  }

  double Ea(double E) {
    this->set_energy(E);

    double Ea_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Ea_ += comp.atoms_bcm * micro_xs.absorption;
    }

    return Ea_;
  }

  double Es(double E) {
    this->set_energy(E);

    double Ea_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Ea_ +=
          comp.atoms_bcm * std::max(micro_xs.total - micro_xs.absorption, 0.);
    }

    return Ea_;
  }

  double Ef(double E) {
    this->set_energy(E);

    double Ef_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Ef_ += comp.atoms_bcm * micro_xs.fission;
    }

    return Ef_;
  }

  double vEf(double E) {
    this->set_energy(E);

    double vEf_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      vEf_ += comp.atoms_bcm * micro_xs.nu_total * micro_xs.fission;
    }

    return vEf_;
  }

  double Eelastic(double E) {
    this->set_energy(E);

    double Eela_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Eela_ += comp.atoms_bcm * micro_xs.elastic;
    }

    return Eela_;
  }

  double Emt(uint32_t mt, double E) {
    this->set_energy(E);

    double Emt_ = 0.;

    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      const auto& eindex = micro_xs.energy_index;
      Emt_ += comp.atoms_bcm * comp.nuclide->reaction_xs(mt, E, eindex);
    }

    return Emt_;
  }

  double heating(double E) {
    this->set_energy(E);

    double Eh_ = 0.;
    // Go through all components in the material
    for (const auto& comp : mat->components()) {
      const auto& micro_xs = this->get_micro_xs(comp.nuclide.get());
      Eh_ += comp.atoms_bcm * micro_xs.heating;
    }
    return Eh_;
  }

  /*
   * This method is used to sample a nuclide based on the traditional branching
   * collision method, where the probability of sampling nuclide i is taken to
   * be $N_i \sigma_i(E) / \Sigma_t(E)$. It can also be used for when performing
   * branchless collisions on the isotope.
   */
  std::pair<const Nuclide*, MicroXSs> sample_nuclide(double E, RNG& rng) {
    this->set_energy(E);

    // First, get total xs
    const double Et_ = this->Et(E);
    const double invs_Et = 1. / Et_;
    const double Ew_ = this->Ew(E);
    const double N_nuclides = static_cast<double>(mat->components().size());

    // Get our random variable
    const double xi = rng();

    // Iterate through all nuclides, untill we find the right one
    double prob_sum = 0.;
    for (const auto& comp : mat->components()) {
      MicroXSs micro_xs = this->get_micro_xs(comp.nuclide.get());
      micro_xs.concentration = comp.atoms_bcm;

      // If we are a noise neutron, we then need to get the copy xs for this
      // nuclide, and adjust the total xs accordingly.
      if (noise_params_) {
        micro_xs.noise_copy = Ew_ / (micro_xs.concentration * N_nuclides);
        micro_xs.total += micro_xs.noise_copy;
      }

      // Get the probability for the nuclide, and add to the cumulative
      const double nuc_prob = invs_Et * micro_xs.concentration * micro_xs.total;
      prob_sum += nuc_prob;

      if (xi <= prob_sum) {
        return {comp.nuclide.get(), micro_xs};
      }
    }

    // We should never get here, but if we do, we just return the last nuclide
    const auto& comp = mat->components().back();
    MicroXSs micro_xs = this->get_micro_xs(comp.nuclide.get());
    micro_xs.concentration = comp.atoms_bcm;

    if (noise_params_) {
      micro_xs.noise_copy = Ew_ / (micro_xs.concentration * N_nuclides);
      micro_xs.total += micro_xs.noise_copy;
    }

    return {comp.nuclide.get(), micro_xs};
  }

  std::pair<const Nuclide*, MicroXSs> sample_branchless_nuclide(
      double E, RNG& rng, BranchlessReaction reaction) {
    this->set_energy(E);

    // First, get the total probability, depending on reaction type
    double sum = 0.;
    for (const auto& comp : mat->components()) {
      const MicroXSs& micro_xs = this->get_micro_xs(comp.nuclide.get());

      switch (reaction) {
        case BranchlessReaction::SCATTER:
          sum += comp.atoms_bcm * (micro_xs.elastic + micro_xs.inelastic);
          break;

        case BranchlessReaction::FISSION:
          sum += comp.atoms_bcm * micro_xs.nu_total * micro_xs.fission;
          break;
      }
    }

    // Get our random variable
    const double xi = rng() * sum;

    // Iterate through all nuclides, untill we find the right one
    double prob_sum = 0.;
    for (const auto& comp : mat->components()) {
      MicroXSs micro_xs = this->get_micro_xs(comp.nuclide.get());
      micro_xs.concentration = comp.atoms_bcm;

      switch (reaction) {
        case BranchlessReaction::SCATTER:
          prob_sum +=
              micro_xs.concentration * (micro_xs.elastic + micro_xs.inelastic);
          break;

        case BranchlessReaction::FISSION:
          prob_sum +=
              micro_xs.concentration * micro_xs.nu_total * micro_xs.fission;
          break;
      }

      if (xi <= prob_sum) {
        return {comp.nuclide.get(), micro_xs};
      }
    }

    // We should never get here, but if we do, we just return the last nuclide
    const auto& comp = mat->components().back();
    MicroXSs micro_xs = this->get_micro_xs(comp.nuclide.get());
    micro_xs.concentration = comp.atoms_bcm;
    return {comp.nuclide.get(), micro_xs};
  }

 private:
  std::optional<NoiseParameters> noise_params_;
  Material* mat;
  double E_;
  boost::unordered_flat_map<const Nuclide*, std::optional<MicroXSs>> xs_;
  boost::unordered_flat_map<uint32_t, std::optional<double>> zaid_to_urr_rand_;

  void clear_xs() {
    // This sets all the MicroXSs objects to nullopt. This typically only
    // called when we have changed energy.
    for (auto& xs : xs_) xs.second = std::nullopt;
  }

  void set_energy(const double& E) {
    if (E_ != E) {
      E_ = E;
      this->clear_xs();
    }
  }

  const MicroXSs& get_micro_xs(const Nuclide* nuc) {
    auto it = xs_.find(nuc);

    // If for some reason there is no entry for this nuclide, we add it.
    if (it == xs_.end()) {
      xs_[nuc] = std::nullopt;
      it = xs_.find(nuc);
    }

    if (it->second.has_value() == false) {
      // We don't have info on this nuclide at this energy yet.
      std::optional<double> urr_rand = std::nullopt;
      const uint32_t zaid = nuc->zaid();

      // If this zaid has URR info, we should use it
      if (settings::use_urr_ptables &&
          zaid_to_urr_rand_.find(zaid) != zaid_to_urr_rand_.end()) {
        urr_rand = zaid_to_urr_rand_[zaid];
      }

      // Get the micro xs and set it
      it->second = nuc->get_micro_xs(E_, urr_rand);
    }

    // Return the micro xs
    return it->second.value();
  }

  const Material::Component& comp(std::size_t i) {
    return mat->components()[i];
  }

  // Function to get the speed of a particle in cm/s
  double speed(double E, std::size_t i) {
    // If we are in CE, it doesn't matter which nuclide we use
    // to get the speed. If we are in MG, then there should
    // only be 1 nuclide, so we always use component 0.
    return mat->components()[0].nuclide->speed(E, i);
  }
};

#endif
