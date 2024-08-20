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
#include <materials/legendre_distribution.hpp>
#include <materials/mg_nuclide.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <cstdint>
#include <sstream>

MGNuclide::MGNuclide(const std::vector<double>& speeds,
                     const std::vector<double>& Et,
                     const std::vector<double>& Ea,
                     const std::vector<double>& Ef,
                     const std::vector<double>& nu_prmpt,
                     const std::vector<double>& nu_dlyd,
                     const std::vector<std::vector<double>>& chi,
                     const std::vector<std::vector<double>>& Es,
                     const std::vector<std::vector<double>>& yield,
                     const std::vector<std::vector<MGAngleDistribution>>& angle,
                     const std::vector<double>& P_dlyd_grp,
                     const std::vector<double>& decay_cnsts,
                     const std::vector<std::vector<double>>& delayed_chi)
    : group_speeds_(speeds),
      Et_(Et),
      Ea_(Ea),
      Ef_(Ef),
      Es_(Ef),  // Initialize with Ef so we have the right size
      nu_prmpt_(nu_prmpt),
      nu_delyd_(nu_dlyd),
      chi_(chi),
      delayed_chi_(delayed_chi),
      Ps_(Es),
      mult_(yield),
      angle_dists_(angle),
      P_delayed_group(P_dlyd_grp),
      delayed_group_decay_constants(decay_cnsts),
      fissile_(false) {
  make_scatter_xs();
  normalize_chi();

  check_sizes();
  check_xs();
  check_fission_data();
  check_dealyed_data();
  check_fissile();
}

void MGNuclide::normalize_chi() {
  // Make sure chi_ is properly sized
  if (chi_.size() != settings::ngroups) {
    fatal_error("chi_.size() is not equal to settings::ngroups.");
  }

  for (std::size_t i = 0; i < settings::ngroups; i++) {
    if (chi_[i].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "chi_[" << i << "].size() is not equal to settings::ngroups.";
      fatal_error(mssg.str());
    }
  }

  for (std::size_t i = 0; i < settings::ngroups; i++) {
    double chi_i = 0.;
    for (std::size_t o = 0; o < settings::ngroups; o++) chi_i += chi_[i][o];
    for (std::size_t o = 0; o < settings::ngroups; o++) chi_[i][o] /= chi_i;
  }

  // Make sure delayed_chi_ is properly sized
  if (delayed_chi_.size() != this->num_delayed_groups()) {
    fatal_error(
        "delayed_chi_.size() is not equal to number of delayed groups.");
  }

  for (std::size_t i = 0; i < this->num_delayed_groups(); i++) {
    if (delayed_chi_[i].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "delayed_chi_[" << i
           << "].size() is not equal to settings::ngroups.";
      fatal_error(mssg.str());
    }
  }

  for (std::size_t i = 0; i < this->num_delayed_groups(); i++) {
    double chi_i = 0.;
    for (std::size_t o = 0; o < settings::ngroups; o++)
      chi_i += delayed_chi_[i][o];
    for (std::size_t o = 0; o < settings::ngroups; o++)
      delayed_chi_[i][o] /= chi_i;
  }
}

void MGNuclide::make_scatter_xs() {
  // Make sure sizes are OK before calculating scatter xs and such
  if (Es_.size() != settings::ngroups) {
    fatal_error("Es_.size() is not equal to settings::ngroups.");
  }

  if (Ps_.size() != settings::ngroups) {
    fatal_error("Ps_.size() is not equal to settings::ngroups.");
  }

  for (std::size_t i = 0; i < settings::ngroups; i++) {
    if (Ps_[i].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Ps_[" << i << "].size() is not equal to settings::ngroups.";
      fatal_error(mssg.str());
    }
  }

  // Calculate total scattering xs.
  for (std::size_t i = 0; i < settings::ngroups; i++) {
    Es_[i] = 0.;
    for (std::size_t o = 0; o < settings::ngroups; o++) Es_[i] += Ps_[i][o];
    for (std::size_t o = 0; o < settings::ngroups; o++) Ps_[i][o] /= Es_[i];
  }
}

void MGNuclide::check_sizes() const {
  if (Et_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of Et_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (Ea_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of Ea_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (Ef_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of Ef_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (Es_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of Es_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (nu_prmpt_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of nu_prmpt_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (nu_delyd_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of nu_delyd_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  if (chi_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of chi_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (chi_[ei].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Size of chi_[" << ei << "] is not " << settings::ngroups;
      mssg << " in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }
  }

  if (Ps_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of Ps_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (Ps_[ei].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Size of Ps_[" << ei << "] is not " << settings::ngroups;
      mssg << " in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }
  }

  if (mult_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of mult_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (mult_[ei].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Size of mult_[" << ei << "] is not " << settings::ngroups;
      mssg << " in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }
  }

  if (angle_dists_.size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Size of angle_dists_ is not " << settings::ngroups;
    mssg << " in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }

  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (angle_dists_[ei].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Size of angle_dists_[" << ei << "] is not " << settings::ngroups;
      mssg << " in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }
  }

  if (P_delayed_group.size() != delayed_group_decay_constants.size()) {
    std::stringstream mssg;
    mssg << "Size of P_delayed_group does not matche the size of ";
    mssg << "delayed_group_decay_constants in nuclide " << this->id() << ".";
    fatal_error(mssg.str());
  }
}

void MGNuclide::check_xs() const {
  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
      // Make sure scattering xs component is positive
      if (Ps_[ei][eo] < 0.) {
        std::stringstream mssg;
        mssg << "Ps_[" << ei << "][" << eo << "] in nuclide " << this->id();
        mssg << " is negative.";
        fatal_error(mssg.str());
      }

      // Make sure yield is positive
      if (mult_[ei][eo] <= 0.) {
        std::stringstream mssg;
        mssg << "mult_[" << ei << "][" << eo << "] in nuclide " << this->id();
        mssg << " is <= 0.";
        fatal_error(mssg.str());
      }
    }

    if (Ea_[ei] < 0.) {
      std::stringstream mssg;
      mssg << "Ea_[" << ei << "] is negative in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }

    if (Ef_[ei] < 0.) {
      std::stringstream mssg;
      mssg << "Ef_[" << ei << "] is negative in nuclide " << this->id() << ".";
      fatal_error(mssg.str());
    }

    if (Ef_[ei] > Ea_[ei]) {
      std::stringstream mssg;
      mssg << "Ef_[" << ei << "] > Ea_[" << ei << "] in nuclide " << this->id()
           << ".";
      fatal_error(mssg.str());
    }

    double diff = Et_[ei] - (Es_[ei] + Ea_[ei]);
    if (std::abs(diff) / Et_[ei] > 0.001) {
      std::stringstream mssg;
      mssg << "In nuclide " << this->id() << ", Es + Ea != Et ";
      mssg << " for group " << ei << ".";
      fatal_error(mssg.str());
    }
  }
}

void MGNuclide::check_fission_data() const {
  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (nu_prmpt_[ei] < 0.) {
      std::stringstream mssg;
      mssg << "nu_prmpt_[" << ei << "] is negative in nuclide ";
      mssg << this->id() << ".";
      fatal_error(mssg.str());
    }

    if (nu_delyd_[ei] < 0.) {
      std::stringstream mssg;
      mssg << "nu_delyd_[" << ei << "] is negative in nuclide " << this->id()
           << ".";
      fatal_error(mssg.str());
    }

    for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
      if (chi_[ei][eo] < 0.) {
        std::stringstream mssg;
        mssg << "chi_[" << ei << "][" << eo << "] is negative in nuclide ";
        mssg << this->id() << ".";
        fatal_error(mssg.str());
      }
    }
  }
}

void MGNuclide::check_dealyed_data() const {
  for (std::size_t i = 0; i < P_delayed_group.size(); i++) {
    if (P_delayed_group[i] < 0.) {
      std::stringstream mssg;
      mssg << "P_delayed_group[" << i << "] is negative in nuclide ";
      mssg << this->id() << ".";
      fatal_error(mssg.str());
    }

    if (delayed_group_decay_constants[i] < 0.) {
      std::stringstream mssg;
      mssg << "delayed_group_decay_constants[" << i
           << "] is negative in nuclide ";
      mssg << this->id() << ".";
      fatal_error(mssg.str());
    }
  }
}

void MGNuclide::check_fissile() {
  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    if (Ef_[ei] != 0. && (nu_prmpt_[ei] + nu_delyd_[ei]) > 0.) {
      fissile_ = true;
      break;
    }
  }
}

bool MGNuclide::fissile() const { return fissile_; }

bool MGNuclide::has_urr() const { return false; }

double MGNuclide::total_xs(double /*E_in*/, std::size_t i) const {
  return Et_[i];
}

double MGNuclide::disappearance_xs(double /*E_in*/, std::size_t i) const {
  return Ea_[i] - Ef_[i];
}

double MGNuclide::fission_xs(double /*E_in*/, std::size_t i) const {
  return Ef_[i];
}

double MGNuclide::nu_total(double /*E_in*/, std::size_t i) const {
  double nu_tot = nu_prmpt_[i];

  if (!nu_delyd_.empty()) nu_tot += nu_delyd_[i];

  return nu_tot;
}

double MGNuclide::nu_prompt(double /*E_in*/, std::size_t i) const {
  return nu_prmpt_[i];
}

double MGNuclide::nu_delayed(double /*E_in*/, std::size_t i) const {
  if (nu_delyd_.empty()) return 0.;

  return nu_delyd_[i];
}

double MGNuclide::reaction_xs(uint32_t /*mt*/, double /*E_in*/,
                              std::size_t /*i*/) const {
  return 0.;
}

double MGNuclide::elastic_xs(double /*E_in*/, std::size_t i) const {
  return Es_[i];
}

double MGNuclide::heating(double /*E_in*/, std::size_t /*i*/) const {
  return 0.;   
}

std::size_t MGNuclide::energy_grid_index(double E) const {
  std::size_t i = 0;

  for (i = 0; i < settings::energy_bounds.size() - 1; i++) {
    if (settings::energy_bounds[i] <= E && E < settings::energy_bounds[i + 1]) {
      break;
    }
  }

  return i;
}

MicroXSs MGNuclide::get_micro_xs(double E,
                                 std::optional<double> /*urr_rand*/) const {
  MicroXSs xs;
  xs.energy = E;
  xs.energy_index = this->energy_grid_index(E);

  xs.total = this->total_xs(E, xs.energy_index);
  xs.fission = this->fission_xs(E, xs.energy_index);
  xs.absorption = xs.fission + this->disappearance_xs(E, xs.energy_index);
  xs.elastic = this->elastic_xs(E, xs.energy_index);
  xs.inelastic = 0.;  // MG has no inelastic
  xs.nu_total = this->nu_total(E, xs.energy_index);
  xs.nu_delayed = this->nu_delayed(E, xs.energy_index);
  xs.concentration = 0.;  // We set this as zero for now
  xs.noise_copy = 0.;     // We also leave this as zero
  xs.heating = this->heating(E, xs.energy_index);
  return xs;
}

std::size_t MGNuclide::num_delayed_groups() const {
  return P_delayed_group.size();
}

double MGNuclide::delayed_group_constant(std::size_t g) const {
  return delayed_group_decay_constants[g];
}

double MGNuclide::delayed_group_probability(std::size_t g, double /*E*/) const {
  return P_delayed_group[g];
}

double MGNuclide::max_energy() const { return settings::energy_bounds.back(); }

double MGNuclide::min_energy() const { return settings::energy_bounds.front(); }

double MGNuclide::speed(double /*E*/, std::size_t i) const {
  return group_speeds_[i];
}

uint32_t MGNuclide::zaid() const { return this->id(); }

double MGNuclide::awr() const {
  std::stringstream mssg;
  mssg << "MGNuclide::awr should never be called. Something bad has happened.";
  fatal_error(mssg.str());
  return 1.;
}

ScatterInfo MGNuclide::sample_scatter(double /*Ein*/, const Direction& u,
                                      const MicroXSs& micro_xs,
                                      RNG& rng) const {
  // Change particle energy
  std::size_t ei = rng.discrete(Ps_[micro_xs.energy_index]);
  double E_out =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Change direction
  std::pair<double, double> mu_and_weight_modifier =
      angle_dists_[micro_xs.energy_index][ei].sample_mu(rng);
  double mu = mu_and_weight_modifier.first;

  double phi = 2. * PI * rng();
  Direction u_out = rotate_direction(u, mu, phi);

  ScatterInfo info;
  info.energy = E_out;
  info.direction = u_out;
  info.mt = 2;
  info.weight_modifier = mu_and_weight_modifier.second;
  return info;
}

ScatterInfo MGNuclide::sample_scatter_mt(uint32_t /*mt*/, double Ein,
                                         const Direction& u, std::size_t i,
                                         RNG& rng) const {
  MicroXSs micro_xs;
  micro_xs.energy_index = i;
  return this->sample_scatter(Ein, u, micro_xs, rng);
}

FissionInfo MGNuclide::sample_prompt_fission(double /*Ein*/, const Direction& u,
                                             std::size_t i, RNG& rng) const {
  // First we sample the the energy index
  std::size_t ei = rng.discrete(chi_[i]);

  // Put fission energy in middle of sampled bin
  double E_out =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * rng() - 1.;
  double phi = 2. * PI * rng();
  Direction u_out = rotate_direction(u, mu, phi);

  FissionInfo info;
  info.energy = E_out;
  info.direction = u_out;
  info.delayed = false;
  info.precursor_decay_constant = 0.;
  info.delayed_family = 0;
  return info;
}

FissionInfo MGNuclide::sample_delayed_fission(double /*Ein*/,
                                              const Direction& u, std::size_t g,
                                              RNG& rng) const {
  FissionInfo info;

  // Sample fission energy
  std::size_t ei = rng.discrete(delayed_chi_[g]);
  // Put fission energy in middle of sampled bin
  info.energy =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Sample direction
  double mu = 2. * rng() - 1.;
  double phi = 2. * PI * rng();
  info.direction = rotate_direction(u, mu, phi);

  // Set delayed group info
  info.delayed = true;
  info.delayed_family = static_cast<uint32_t>(g);
  info.precursor_decay_constant = this->delayed_group_decay_constants[g];

  return info;
}

FissionInfo MGNuclide::sample_fission(double /*Ein*/, const Direction& u,
                                      std::size_t energy_index, double Pdelayed,
                                      RNG& rng) const {
  FissionInfo info;

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * rng() - 1.;
  double phi = 2. * PI * rng();
  Direction uout = rotate_direction(u, mu, phi);

  info.direction = uout;
  info.delayed = false;
  info.delayed_family = 0;
  info.precursor_decay_constant = 0.;

  // Variable to hold outgoing energy index
  std::size_t ei;

  // Next, we need to see if this is a delayed neutron or not.
  if (rng() < Pdelayed) {
    // We have a delayed neutron. We now need to select a delayed group
    // and get the group decay constant.
    std::size_t dgrp = rng.discrete(P_delayed_group);
    double lambda = delayed_group_decay_constants[dgrp];

    ei = rng.discrete(delayed_chi_[dgrp]);

    info.delayed = true;
    info.delayed_family = static_cast<uint32_t>(dgrp);
    info.precursor_decay_constant = lambda;
  } else {
    ei = rng.discrete(chi_[energy_index]);
  }

  // Put fission energy in middle of sampled bin
  info.energy =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  return info;
}

void get_legendre_moment(
    const YAML::Node& mat, uint32_t id, std::size_t l,
    std::vector<std::vector<LegendreDistribution>>& angles) {
  std::string key = "P" + std::to_string(l);

  if (mat[key]) {
    if (!mat[key].IsSequence() || mat[key].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Invalid " << key << " matrix entry in material " << id << ".";
      fatal_error(mssg.str());
    }

    // Go through all rows (incoming energy groups)
    for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
      // Make sure the row is the right size
      if (!mat[key][ei].IsSequence() ||
          mat[key][ei].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Row " << ei << " of the " << key << " matrix for material "
             << id;
        mssg << " is invalid.";
        fatal_error(mssg.str());
      }

      // Go through all columns (outgoing energy groups)
      for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
        double coeff = mat[key][ei][eo].as<double>();
        angles[ei][eo].set_moment(l, coeff);
      }
    }
  }
}

std::shared_ptr<MGNuclide> make_mg_nuclide(const YAML::Node& mat, uint32_t id) {
  //===========================================================================
  // Get Total XS
  std::vector<double> Et;
  if (!mat["total"] || !mat["total"].IsSequence() ||
      mat["total"].size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Invalid total xs entry in material " << id << ".";
    fatal_error(mssg.str());
  }
  Et = mat["total"].as<std::vector<double>>();

  //===========================================================================
  // Get Absorption XS
  std::vector<double> Ea;
  if (!mat["absorption"] || !mat["absorption"].IsSequence() ||
      mat["absorption"].size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Invalid absorption xs entry in material " << id << ".";
    fatal_error(mssg.str());
  }
  Ea = mat["absorption"].as<std::vector<double>>();

  //===========================================================================
  // Get the scatter matrix
  std::vector<std::vector<double>> Es(
      settings::ngroups, std::vector<double>(settings::ngroups, 0.));
  if (!mat["scatter"] || !mat["scatter"].IsSequence() ||
      mat["scatter"].size() != settings::ngroups) {
    std::stringstream mssg;
    mssg << "Invalid scatter matrix entry in material " << id << ".";
    fatal_error(mssg.str());
  }

  // Go through all rows (incoming energy groups)
  for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
    // Make sure the row is the right size
    if (!mat["scatter"][ei].IsSequence() ||
        mat["scatter"][ei].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Row " << ei << " of the scatter matrix for material " << id;
      mssg << " is invalid.";
      fatal_error(mssg.str());
    }

    // Go through all columns (outgoing energy groups)
    for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
      Es[ei][eo] = mat["scatter"][ei][eo].as<double>();

      if (Es[ei][eo] < 0.) {
        std::stringstream mssg;
        mssg << "Negative scattering component at row " << ei;
        mssg << ", column " << eo << " of the scattering matrix";
        mssg << " for material " << id << ".";
        fatal_error(mssg.str());
      }
    }

    // Check the ingroup scattering. If any ingroup scattering is zero,
    // we need can't use virtual collisions in exact cancellation.
    if (Es[ei][ei] == 0.) settings::use_virtual_collisions = false;
  }

  //===========================================================================
  // Get Scattering yields (If Present)
  std::vector<std::vector<double>> yields(
      settings::ngroups, std::vector<double>(settings::ngroups, 1.));
  if (mat["yields"]) {
    if (!mat["yields"].IsSequence() ||
        mat["yields"].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Invalid yields matrix entry in material " << id << ".";
      fatal_error(mssg.str());
    }

    // Go through all rows (incoming energy groups)
    for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
      // Make sure the row is the right size
      if (!mat["yields"][ei].IsSequence() ||
          mat["yields"][ei].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Row " << ei << " of the yields matrix for material " << id;
        mssg << " is invalid.";
        fatal_error(mssg.str());
      }

      // Go through all columns (outgoing energy groups)
      for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
        yields[ei][eo] = mat["yields"][ei][eo].as<double>();

        if (yields[ei][eo] < 0.) {
          std::stringstream mssg;
          mssg << "Negative scattering yield at row " << ei;
          mssg << ", column " << eo << " of the scattering matrix";
          mssg << " for material " << id << ".";
          fatal_error(mssg.str());
        }
      }
    }
  }

  //===========================================================================
  // Get Angular Distributions (If Present)
  std::vector<std::vector<LegendreDistribution>> legendre_angles(
      settings::ngroups, std::vector<LegendreDistribution>(settings::ngroups));

  // Get up to 5 legendre moments, if they are present.
  for (std::size_t l = 1; l <= 5; l++) {
    get_legendre_moment(mat, id, l, legendre_angles);
  }

  // Make sure all distributions are positive everywhere.
  for (std::size_t i = 0; i < settings::ngroups; i++) {
    for (std::size_t o = 0; o < settings::ngroups; o++) {
      if (legendre_angles[i][o].positive_over_domain() == false) {
        std::stringstream mssg;
        mssg << "Angle distribution for scattering from group " << i;
        mssg << " to group " << o << " for material " << id;
        mssg << " has negative regions.";
        warning(mssg.str());
      }
    }
  }

  std::vector<std::vector<MGAngleDistribution>> angles(
      settings::ngroups, std::vector<MGAngleDistribution>(settings::ngroups));

  for (std::size_t i = 0; i < settings::ngroups; i++) {
    for (std::size_t o = 0; o < settings::ngroups; o++) {
      angles[i][o] = legendre_angles[i][o].linearize();
    }

    // Check the in-scattering distributions. If zero for mu = 1, then we
    // can't use virtual collisions in exact cancellation.
    if (angles[i][i].pdf(1.) == 0.) settings::use_virtual_collisions = false;
  }

  //===========================================================================
  // Get Fission XS if Present
  std::vector<double> Ef(settings::ngroups, 0.);
  bool fissile = false;
  if (mat["fission"]) {
    if (!mat["fission"].IsSequence() ||
        mat["fission"].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Invalid fission xs entry in material " << id << ".";
      fatal_error(mssg.str());
    }
    Ef = mat["fission"].as<std::vector<double>>();

    for (const auto& xs_f : Ef) {
      if (xs_f < 0.) {
        std::stringstream mssg;
        mssg << "Negative fission xs in material " << id << ".";
        fatal_error(mssg.str());
      }

      if (xs_f > 0.) fissile = true;
    }
  }

  //===========================================================================
  // Get Nu if we have fission
  std::vector<double> nu_prmpt(settings::ngroups, 0.);
  std::vector<double> nu_dlyd(settings::ngroups, 0.);
  if (fissile) {
    if (mat["nu"]) {
      // Treat nu as nu_prmpt
      if (!mat["nu"].IsSequence() || mat["nu"].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Invalid nu entry in material " << id << ".";
        fatal_error(mssg.str());
      }
      nu_prmpt = mat["nu"].as<std::vector<double>>();
    } else if (mat["nu_prompt"] && mat["nu_delayed"]) {
      // Prompt
      if (!mat["nu_prompt"].IsSequence() ||
          mat["nu_prompt"].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Invalid nu_prompt entry in material " << id << ".";
        fatal_error(mssg.str());
      }
      nu_prmpt = mat["nu_prompt"].as<std::vector<double>>();

      // Delayed
      if (!mat["nu_delayed"].IsSequence() ||
          mat["nu_delayed"].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Invalid nu_delayed entry in material " << id << ".";
        fatal_error(mssg.str());
      }
      nu_dlyd = mat["nu_delayed"].as<std::vector<double>>();
    } else {
      if (mat["nu_prompt"]) {
        // No nu_delayed data is given. This is bad
        std::stringstream mssg;
        mssg << "No nu_delayed data is provided in material " << id << ".";
        fatal_error(mssg.str());
      } else if (mat["nu_delayed"]) {
        // No nu_prompt data is given. This is bad
        std::stringstream mssg;
        mssg << "No nu_prompt data is provided in material " << id << ".";
        fatal_error(mssg.str());
      } else {
        // No nu data is given at all. This is really bad
        std::stringstream mssg;
        mssg << "No nu data is provided in material " << id << ".";
        fatal_error(mssg.str());
      }
    }
  }

  //===========================================================================
  // Get Chi Matrix if we have fission
  std::vector<std::vector<double>> chi(
      settings::ngroups, std::vector<double>(settings::ngroups, 0.));
  if (fissile) {
    if (!mat["chi"] || !mat["chi"].IsSequence()) {
      std::stringstream mssg;
      mssg << "Invalid chi matrix entry in material " << id << ".";
      fatal_error(mssg.str());
    }

    if (mat["chi"].size() != 1 && mat["chi"].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Invalid chi matrix entry in material " << id << ".";
      fatal_error(mssg.str());
    }

    const bool chi_matrix = (mat["chi"].size() == settings::ngroups);
    // Set flag needed for ExactMGCancelator
    if (chi_matrix) settings::chi_matrix = true;

    if (chi_matrix) {
      // Make sure all rows have the right length
      for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
        if (!mat["chi"][ei].IsSequence() ||
            mat["chi"][ei].size() != settings::ngroups) {
          std::stringstream mssg;
          mssg << "Invalid length for chi[" << ei << "] in material " << id
               << ".";
          fatal_error(mssg.str());
        }

        for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
          chi[ei][eo] = mat["chi"][ei][eo].as<double>();

          if (chi[ei][eo] < 0.) {
            std::stringstream mssg;
            mssg << "chi[" << ei << "][" << eo << "] is negative in material "
                 << id << ".";
            fatal_error(mssg.str());
          }
        }
      }
    } else {
      if (!mat["chi"][0].IsSequence() ||
          mat["chi"][0].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "Invalid length for chi[0] in material " << id << ".";
        fatal_error(mssg.str());
      }

      for (std::size_t ei = 0; ei < settings::ngroups; ei++) {
        for (std::size_t eo = 0; eo < settings::ngroups; eo++) {
          chi[ei][eo] = mat["chi"][0][eo].as<double>();

          if (chi[ei][eo] < 0.) {
            std::stringstream mssg;
            mssg << "chi[" << ei << "][" << eo << "] is negative in material "
                 << id << ".";
            fatal_error(mssg.str());
          }
        }
      }
    }
  }

  //===========================================================================
  // Get the delayed group info
  std::vector<double> P_delayed_grp;
  std::vector<double> delayed_constants;
  std::vector<std::vector<double>> delayed_chi;
  if (mat["delayed_groups"] && mat["delayed_groups"].IsMap()) {
    // Get the probability of each group
    if (!mat["delayed_groups"]["probabilities"] ||
        !mat["delayed_groups"]["probabilities"].IsSequence()) {
      std::stringstream mssg;
      mssg << "No probabilities entry in delayed_groups for material " << id
           << ".";
      fatal_error(mssg.str());
    }
    P_delayed_grp =
        mat["delayed_groups"]["probabilities"].as<std::vector<double>>();

    // Get the decay constant of each group
    if (!mat["delayed_groups"]["constants"] ||
        !mat["delayed_groups"]["constants"].IsSequence()) {
      std::stringstream mssg;
      mssg << "No constants entry in delayed_groups for material " << id << ".";
      fatal_error(mssg.str());
    }
    delayed_constants =
        mat["delayed_groups"]["constants"].as<std::vector<double>>();

    // Make sure both have the same size
    if (P_delayed_grp.size() != delayed_constants.size()) {
      std::stringstream mssg;
      mssg << "In delayed_groups entry for material " << id << ", ";
      mssg << "probabilities and constants entries have different sizes.";
      fatal_error(mssg.str());
    }

    // Get chi for delayed groups
    if (!mat["delayed_groups"]["chi"] ||
        mat["delayed_groups"]["chi"].IsSequence() == false) {
      std::stringstream mssg;
      mssg << "No chi entry in delayed groups for material " << id << ".";
      fatal_error(mssg.str());
    }
    delayed_chi =
        mat["delayed_groups"]["chi"].as<std::vector<std::vector<double>>>();

    // Make sure delayed chi has the right size
    if (delayed_chi.size() != P_delayed_grp.size()) {
      std::stringstream mssg;
      mssg << "In delayed_groups entry for material " << id << ", ";
      mssg << "probabilities and chi entries have different sizes.";
      fatal_error(mssg.str());
    }

    // Make sure each delayed group spectrum has the correct number of entries
    for (std::size_t d = 0; d < delayed_chi.size(); d++) {
      if (delayed_chi[d].size() != settings::ngroups) {
        std::stringstream mssg;
        mssg << "In delayed_groups entry for material " << id << ", ";
        mssg << "the length of the chi entry for delayed family " << d;
        mssg << " does not agree with the number of energy groups.";
        fatal_error(mssg.str());
      }

      // Make sure chi is positive and not all zeros
      double chi_d_sum = 0.;
      for (std::size_t g = 0; g < settings::ngroups; g++) {
        chi_d_sum += delayed_chi[d][g];

        if (delayed_chi[d][g] < 0.) {
          std::stringstream mssg;
          mssg << "In delayed_groups entry for material " << id << ", ";
          mssg << "chi[" << d << "][" << g << "] is less than zero.";
          fatal_error(mssg.str());
        }
      }
      if (chi_d_sum == 0.) {
        std::stringstream mssg;
        mssg << "In delayed_groups entry for material " << id << ", ";
        mssg << "chi for delayed family " << d << " is all zeros.";
        fatal_error(mssg.str());
      }
    }
  } else if (mat["delayed_groups"]) {
    std::stringstream mssg;
    mssg << "Invalid delayed_groups entry in material " << id << ".";
    fatal_error(mssg.str());
  }

  //===========================================================================
  // Get the group speeds
  std::vector<double> grp_speeds(settings::ngroups, 1.);
  if (mat["group-speeds"]) {
    if (!mat["group-speeds"].IsSequence() ||
        mat["group-speeds"].size() != settings::ngroups) {
      std::stringstream mssg;
      mssg << "Invalid group-speeds entry in material " << id << ".";
      fatal_error(mssg.str());
    }
    grp_speeds = mat["group-speeds"].as<std::vector<double>>();
  } else if (settings::sim_mode == settings::SimMode::NOISE) {
    std::stringstream mssg;
    mssg << "Missing group-speeds entry in material " << id << ".";
    fatal_error(mssg.str());
  }

  // We should have all info ! Now we can return the nuclide
  auto nuclide_to_return = std::make_shared<MGNuclide>(
      grp_speeds, Et, Ea, Ef, nu_prmpt, nu_dlyd, chi, Es, yields, angles,
      P_delayed_grp, delayed_constants, delayed_chi);

  // Record nuclide in nuclides map
  nuclides[nuclide_to_return->id()] = nuclide_to_return;

  return nuclide_to_return;
}
