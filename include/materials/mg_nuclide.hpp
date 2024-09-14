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
#ifndef MG_NUCLIDE_H
#define MG_NUCLIDE_H

#include <materials/mg_angle_distribution.hpp>
#include <materials/nuclide.hpp>

#include <yaml-cpp/yaml.h>

#include <cstdint>
#include <vector>

class MGNuclide : public Nuclide {
 public:
  MGNuclide(const std::vector<double>& speeds, const std::vector<double>& Et,
            const std::vector<double>& Ea, const std::vector<double>& Ef,
            const std::vector<double>& nu_prmpt,
            const std::vector<double>& nu_dlyd,
            const std::vector<std::vector<double>>& chi,
            const std::vector<std::vector<double>>& Es,
            const std::vector<std::vector<double>>& yield,
            const std::vector<std::vector<MGAngleDistribution>>& angle,
            const std::vector<double>& P_dlyd_grp,
            const std::vector<double>& decay_cnsts,
            const std::vector<std::vector<double>>& delayed_chi);

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

  // Methods unique to MGNuclide, for exact MG cancellation.
  const std::vector<double>& group_speeds() const { return group_speeds_; }
  const std::vector<double>& Et() const { return Et_; }
  const std::vector<double>& Ea() const { return Ea_; }
  const std::vector<double>& Ef() const { return Ef_; }
  const std::vector<double>& Es() const { return Es_; }
  const std::vector<double>& nu_prmpt() const { return nu_prmpt_; }
  const std::vector<double>& nu_dlyd() const { return nu_delyd_; }
  const std::vector<std::vector<double>>& chi() const { return chi_; }
  const std::vector<std::vector<double>>& Ps() const { return Ps_; }
  const std::vector<std::vector<double>>& mult() const { return mult_; }
  const std::vector<std::vector<MGAngleDistribution>>& angles() const {
    return angle_dists_;
  }
  const std::vector<double>& P_dlyd_grp() const { return P_delayed_group; }
  const std::vector<double>& dlyd_grp_decay_const() const {
    return delayed_group_decay_constants;
  }

 private:
  std::vector<double> group_speeds_;
  std::vector<double> Et_;
  std::vector<double> Ea_;
  std::vector<double> Ef_;
  std::vector<double> Es_;
  std::vector<double> nu_prmpt_;
  std::vector<double> nu_delyd_;
  std::vector<std::vector<double>> chi_;
  std::vector<std::vector<double>> delayed_chi_;
  std::vector<std::vector<double>> Ps_;
  std::vector<std::vector<double>> mult_;
  std::vector<std::vector<MGAngleDistribution>> angle_dists_;
  std::vector<double> P_delayed_group;
  std::vector<double> delayed_group_decay_constants;
  bool fissile_ = false;

  void make_scatter_xs();
  void normalize_chi();
  void check_sizes() const;
  void check_xs() const;
  void check_fission_data() const;
  void check_dealyed_data() const;
  void check_fissile();
};

std::shared_ptr<MGNuclide> make_mg_nuclide(const YAML::Node& mat, uint32_t id);

#endif
