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
#include <utils/constants.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <memory>
#include <random>

namespace settings {

SimMode sim_mode = SimMode::KEFF;

uint32_t ngroups = 0;

bool plotting_mode = false;

Timer alpha_omega_timer;
double max_time = INF;

double min_energy = 0.;       // [MeV]
double max_energy = 100000.;  // [MeV]
double target_at_rest_threshold = 400.;
bool use_urr_ptables = true;

EnergyMode energy_mode = EnergyMode::CE;

uint64_t rng_seed = 19073486328125;
uint64_t rng_stride = 152917;
RNG rng;
bool rng_stride_warnings = false;

// double wgt_cutoff = 0.8;  // Needs to be high for noise !
double wgt_cutoff = 0.25;
double wgt_survival = 1.0;
double wgt_split = 2.0;

// Energy bounds for multi-group mode
std::vector<double> energy_bounds{};
bool chi_matrix = false;
bool use_virtual_collisions = true;
std::size_t group(double E) {
  if (energy_bounds.size() <= 1) return 0;

  for (std::size_t g = 0; g < energy_bounds.size() - 1; g++) {
    if (energy_bounds[g] <= E && E <= energy_bounds[g + 1]) return g;
  }

  // Should never get here
  return 0;
}

std::unique_ptr<NDDirectory> nd_directory = nullptr;
std::string nd_directory_fname;
bool use_dbrc = true;
std::vector<std::string> dbrc_nuclides;
TempInterpolation temp_interpolation = TempInterpolation::Linear;

void initialize_global_rng() {
  rng.seed(rng_seed);
  rng.set_stream(2);
}

void write_settings_to_output() {
  if (mpi::rank != 0) return;

  // Get reference to HDF5 file
  H5::File& h5 = Output::instance().h5();

  // Energy Mode
  switch (energy_mode) {
    case EnergyMode::CE:
      h5.createAttribute<std::string>("energy-mode", "continuous-energy");
      break;
    case EnergyMode::MG:
      h5.createAttribute<std::string>("energy-mode", "multi-group");
      break;
  }

  // MG CT specific
  if (energy_mode == EnergyMode::MG) {
    h5.createAttribute("ngroups", ngroups);

    h5.createAttribute("energy-bounds", energy_bounds);
  }

  if (energy_mode == EnergyMode::CE) {
    switch (temp_interpolation) {
      case TempInterpolation::Exact:
        h5.createAttribute<std::string>("temperature-interpolation", "exact");
        break;
      case TempInterpolation::Nearest:
        h5.createAttribute<std::string>("temperature-interpolation", "nearest");
        break;
      case TempInterpolation::Linear:
        h5.createAttribute<std::string>("temperature-interpolation", "linear");
        break;
    }

    h5.createAttribute<std::string>("nuclear-data", nd_directory_fname);

    h5.createAttribute<bool>("use-dbrc", use_dbrc);

    if (dbrc_nuclides.size() > 0) {
      h5.createAttribute("dbrc-nuclides", dbrc_nuclides);
    }

    h5.createAttribute<bool>("use-urr-ptables", use_urr_ptables);
  }

  // Common bits
  h5.createAttribute("max-run-time", max_time);

  h5.createAttribute("rng-stride-warnings", rng_stride_warnings);

  h5.createAttribute("seed", rng_seed);

  h5.createAttribute("stride", rng_stride);

  h5.createAttribute("wgt-cutoff", wgt_cutoff);

  h5.createAttribute("wgt-survival", wgt_survival);

  h5.createAttribute("wgt-split", wgt_split);
}
}  // namespace settings
