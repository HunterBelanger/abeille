/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
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
#ifndef CE_NUCLIDE_H
#define CE_NUCLIDE_H

#include <materials/nuclide.hpp>

#include <PapillonNDL/st_neutron.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>

#include <cstdint>
#include <memory>
#include <vector>

class CENuclide : public Nuclide {
 public:
  CENuclide(const std::shared_ptr<pndl::STNeutron>& ce,
            const std::shared_ptr<pndl::STThermalScatteringLaw>& tsl);

  bool fissile() const override final;
  bool has_urr() const override final;
  double total_xs(double E_in, std::size_t i) const override final;
  double disappearance_xs(double E_in, std::size_t i) const override final;
  double fission_xs(double E_in, std::size_t i) const override final;
  double nu_total(double E_in, std::size_t i) const override final;
  double nu_prompt(double E_in, std::size_t i) const override final;
  double nu_delayed(double E_in, std::size_t i) const override final;
  double reaction_xs(uint32_t mt, double E_in, size_t i) const override final;
  double elastic_xs(double E_in, std::size_t i) const override final;
  double heating(double E_in, std::size_t i) const override final;
  std::size_t energy_grid_index(double E) const override final;
  MicroXSs get_micro_xs(double E, std::optional<double> urr_rand =
                                      std::nullopt) const override final;

  std::size_t num_delayed_groups() const override final;
  double delayed_group_constant(std::size_t g) const override final;
  double delayed_group_probability(std::size_t g,
                                   double E) const override final;

  ScatterInfo sample_scatter(double Ein, const Direction& u,
                             const MicroXSs& micro_xs,
                             RNG& rng) const override final;
  ScatterInfo sample_scatter_mt(uint32_t mt, double Ein, const Direction& u,
                                std::size_t i, RNG& rng) const override final;
  FissionInfo sample_fission(double Ein, const Direction& u, std::size_t i,
                             double Pdelayed, RNG& rng) const override final;
  FissionInfo sample_prompt_fission(double Ein, const Direction& u,
                                    std::size_t i,
                                    RNG& rng) const override final;
  FissionInfo sample_delayed_fission(double Ein, const Direction& u,
                                     std::size_t g,
                                     RNG& rng) const override final;

  double max_energy() const override final;
  double min_energy() const override final;
  double speed(double E, std::size_t i) const override final;

  uint32_t zaid() const override final;
  double awr() const override final;

  // CENuclide specific options
  const std::shared_ptr<pndl::STNeutron>& cedata() const { return cedata_; }
  const std::shared_ptr<pndl::STThermalScatteringLaw>& tsl() const {
    return tsl_;
  }
  const std::vector<double>& urr_energy_grid() const {
    return cedata_->urr_ptables().energy();
  }

 private:
  std::shared_ptr<pndl::STNeutron> cedata_;
  std::shared_ptr<pndl::STThermalScatteringLaw> tsl_;

  void elastic_scatter(double Ein, const Direction& uin, double& Eout,
                       Direction& uout, RNG& rng) const;

  void thermal_scatter(double Ein, const Direction& uin, double& Eout,
                       Direction& uout, RNG& rng) const;
};

#endif
