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

#include <complex>
#include <cstdint>
#include <functional>
#include <materials/legendre_distribution.hpp>
#include <materials/mg_nuclide.hpp>
#include <sstream>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

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
                     const std::vector<double>& decay_cnsts)
    : group_speeds_(speeds),
      Et_(Et),
      Ea_(Ea),
      Ef_(Ef),
      Es_(Ef),  // Initialize with Ef so we have the right size
      nu_prmpt_(nu_prmpt),
      nu_delyd_(nu_dlyd),
      chi_(chi),
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
                                      pcg32& rng) const {
  // Change particle energy
  std::size_t ei =
      static_cast<std::size_t>(RNG::discrete(rng, Ps_[micro_xs.energy_index]));
  double E_out =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Change direction
  double mu = angle_dists_[micro_xs.energy_index][ei].sample_mu(rng);
  double phi = 2. * PI * RNG::rand(rng);
  Direction u_out = rotate_direction(u, mu, phi);

  ScatterInfo info;
  info.energy = E_out;
  info.direction = u_out;
  info.mt = 2;
  return info;
}

ScatterInfo MGNuclide::sample_scatter_mt(uint32_t /*mt*/, double Ein,
                                         const Direction& u, std::size_t i,
                                         pcg32& rng) const {
  MicroXSs micro_xs;
  micro_xs.energy_index = i;
  return this->sample_scatter(Ein, u, micro_xs, rng);
}

FissionInfo MGNuclide::sample_prompt_fission(double /*Ein*/, const Direction& u,
                                             std::size_t i, pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));

  // First we sample the the energy index
  std::size_t ei = static_cast<std::size_t>(RNG::discrete(rng, chi_[i]));

  // Put fission energy in middle of sampled bin
  double E_out =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * RNG::rand(rng) - 1.;
  double phi = 2. * PI * RNG::rand(rng);
  Direction u_out = rotate_direction(u, mu, phi);

  FissionInfo info;
  info.energy = E_out;
  info.direction = u_out;
  info.delayed = false;
  info.precursor_decay_constant = 0.;
  info.delayed_family = 0;
  return info;
}

FissionInfo MGNuclide::sample_delayed_fission(double Ein, const Direction& u,
                                              std::size_t g, pcg32& rng) const {
  FissionInfo info = sample_prompt_fission(Ein, u, g, rng);
  info.delayed_family = static_cast<uint32_t>(g);
  info.precursor_decay_constant = this->delayed_group_decay_constants[g];
  return info;
}

FissionInfo MGNuclide::sample_fission(double /*Ein*/, const Direction& u,
                                      std::size_t energy_index, double Pdelayed,
                                      pcg32& rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));

  FissionInfo info;

  // First we sample the the energy index
  std::size_t ei =
      static_cast<std::size_t>(RNG::discrete(rng, chi_[energy_index]));

  // Put fission energy in middle of sampled bin
  double E_out =
      0.5 * (settings::energy_bounds[ei] + settings::energy_bounds[ei + 1]);

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * RNG::rand(rng) - 1.;
  double phi = 2. * PI * RNG::rand(rng);
  Direction uout = rotate_direction(u, mu, phi);

  info.energy = E_out;
  info.direction = uout;
  info.delayed = false;
  info.delayed_family = 0;
  info.precursor_decay_constant = 0.;

  // Next, we need to see if this is a delayed neutron or not.
  if (rngfunc() < Pdelayed) {
    // We have a delayed neutron. We now need to select a delayed group
    // and get the group decay constant.
    std::size_t dgrp =
        static_cast<std::size_t>(RNG::discrete(rng, P_delayed_group));
    double lambda = delayed_group_decay_constants[dgrp];

    info.delayed = true;
    info.delayed_family = static_cast<uint32_t>(dgrp);
    info.precursor_decay_constant = lambda;
  }

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
  } else if (settings::mode == settings::SimulationMode::NOISE) {
    std::stringstream mssg;
    mssg << "Missing group-speeds entry in material " << id << ".";
    fatal_error(mssg.str());
  }

  // We should have all info ! Now we can return the nuclide
  auto nuclide_to_return = std::make_shared<MGNuclide>(
      grp_speeds, Et, Ea, Ef, nu_prmpt, nu_dlyd, chi, Es, yields, angles,
      P_delayed_grp, delayed_constants);

  // Record nuclide in nuclides map
  nuclides[nuclide_to_return->id()] = nuclide_to_return;

  return nuclide_to_return;
}
