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
#ifndef NUCLIDE_H
#define NUCLIDE_H

#include <noise_source/noise_source.hpp>
#include <simulation/particle.hpp>

#include <cstdint>
#include <map>
#include <optional>
#include <unordered_set>

class Nuclide;
extern std::map<uint32_t, std::shared_ptr<Nuclide>> nuclides;
extern std::unordered_set<uint32_t> zaids_with_urr;

struct MicroXSs {
  double total = 0.;
  double fission = 0.;
  double absorption = 0.;
  double elastic = 0.;
  double inelastic =
      0.;  // Not actually inelastic, just all scattering that isn't elastic
  double nu_total = 0.;
  double nu_delayed = 0.;
  double energy = 0.;
  double concentration = 0.;
  double noise_copy = 0.;
  std::size_t energy_index = 0;
  double heating = 0.;
  bool urr = false;
};

struct ScatterInfo {
  double yield = 1.;
  double energy = 0.;
  Direction direction = Direction(1., 0., 0.);
  uint32_t mt = 0;
  bool thermal = false;
  double weight_modifier = 1.0;
};

struct FissionInfo {
  double energy = 0.;
  Direction direction = Direction(1., 0., 0.);
  bool delayed = false;
  double precursor_decay_constant = 0.;
  uint32_t delayed_family = 0;
};

// This is a general nuclide interface, to permit use with
// the same transport and simulation classes, regardless of
// the energy mode being used.
class Nuclide {
 public:
  Nuclide() : id_(id_counter++) {}
  virtual ~Nuclide() = default;

  virtual bool fissile() const = 0;
  virtual bool has_urr() const = 0;
  virtual double total_xs(double E_in, std::size_t i) const = 0;
  virtual double disappearance_xs(double E_in, std::size_t i) const = 0;
  virtual double fission_xs(double E_in, std::size_t i) const = 0;
  virtual double nu_total(double E_in, std::size_t i) const = 0;
  virtual double nu_prompt(double E_in, std::size_t i) const = 0;
  virtual double nu_delayed(double E_in, std::size_t i) const = 0;
  virtual double reaction_xs(uint32_t mt, double E_in, size_t i) const = 0;
  virtual double elastic_xs(double E_in, std::size_t i) const = 0;
  virtual double heating(double E_in, std::size_t i) const = 0;
  virtual std::size_t energy_grid_index(double E) const = 0;
  virtual MicroXSs get_micro_xs(
      double E, std::optional<double> urr_rand = std::nullopt) const = 0;

  virtual std::size_t num_delayed_groups() const = 0;
  virtual double delayed_group_constant(std::size_t g) const = 0;
  virtual double delayed_group_probability(std::size_t g, double E) const = 0;

  virtual ScatterInfo sample_scatter(double Ein, const Direction& u,
                                     const MicroXSs& micro_xs,
                                     RNG& rng) const = 0;
  virtual ScatterInfo sample_scatter_mt(uint32_t mt, double Ein,
                                        const Direction& u, std::size_t i,
                                        RNG& rng) const = 0;
  virtual FissionInfo sample_fission(double Ein, const Direction& u,
                                     std::size_t i, double Pdelayed,
                                     RNG& rng) const = 0;
  virtual FissionInfo sample_prompt_fission(double Ein, const Direction& u,
                                            std::size_t i, RNG& rng) const = 0;
  // Samples an energy and direction from the delayed spectrum of delayed family
  // g.
  virtual FissionInfo sample_delayed_fission(double Ein, const Direction& u,
                                             std::size_t g, RNG& rng) const = 0;

  virtual double max_energy() const = 0;
  virtual double min_energy() const = 0;
  virtual double speed(double E, std::size_t i) const = 0;

  virtual uint32_t zaid() const = 0;
  virtual double awr() const = 0;

  uint32_t id() const { return id_; }

 private:
  static uint32_t id_counter;
  uint32_t id_;
};

#endif
