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
int nparticles = 100000;  // Number of particles per batch / generation
int ngenerations = 120;
int nignored = 20;
int nskip =
    10;  // Number of power iteration generations between each noise batch

uint32_t ngroups = 0;
int n_cancel_noise_gens = INF_INT;

bool plotting_mode = false;

Timer alpha_omega_timer;
double max_time = INF;

double min_energy = 0.;       // [MeV]
double max_energy = 100000.;  // [MeV]

double target_at_rest_threshold = 400.;

SimulationMode mode = SimulationMode::K_EIGENVALUE;
TrackingMode tracking = TrackingMode::SURFACE_TRACKING;
EnergyMode energy_mode = EnergyMode::CE;

uint64_t rng_seed = 19073486328125;
uint64_t rng_stride = 152917;
pcg32 rng;

// double wgt_cutoff = 0.8;  // Needs to be high for noise !
double wgt_cutoff = 0.25;
double wgt_survival = 1.0;
double wgt_split = 2.0;

double w_noise = -1.;
double eta = 1.;   // Used in noise transport
double keff = 1.;  // Used in noise transport

bool use_urr_ptables = true;

bool converged = false;

bool pair_distance_sqrd = false;
bool families = false;
bool empty_entropy_bins = false;

bool regional_cancellation = false;
bool regional_cancellation_noise = false;

bool inner_generations = true;
bool normalize_noise_source = true;
bool rng_stride_warnings = false;

bool load_source_file = false;

bool branchless_splitting = false;
bool branchless_combing = true;
bool branchless_material = true;

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
std::vector<double> sample_xs_ratio{};

std::string output_file_name = "output.txt";
std::string in_source_file_name = "";

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

  // Simulation Mode
  switch (mode) {
    case SimulationMode::K_EIGENVALUE:
      h5.createAttribute<std::string>("mode", "k-eigenvalue");
      break;
    case SimulationMode::BRANCHLESS_K_EIGENVALUE:
      h5.createAttribute<std::string>("mode", "branchless-k-eigenvalue");
      break;
    case SimulationMode::FIXED_SOURCE:
      h5.createAttribute<std::string>("mode", "fixed-source");
      break;
    case SimulationMode::MODIFIED_FIXED_SOURCE:
      h5.createAttribute<std::string>("mode", "modified-fixed-source");
      break;
    case SimulationMode::NOISE:
      h5.createAttribute<std::string>("mode", "noise");
      break;
  }

  // Branchless specific things
  if (mode == SimulationMode::BRANCHLESS_K_EIGENVALUE) {
    h5.createAttribute<bool>("branchless-combing", branchless_combing);
    h5.createAttribute<bool>("branchless-splitting", branchless_splitting);
    h5.createAttribute<bool>("branchless-material", branchless_material);
  }

  // Noise specific things
  if (mode == SimulationMode::NOISE) {
    h5.createAttribute("nskip", nskip);

    h5.createAttribute("noise-angular-frequency", w_noise);

    h5.createAttribute("keff", keff);

    h5.createAttribute<bool>("inner-generations", inner_generations);

    h5.createAttribute<bool>("normalize-noise-source", normalize_noise_source);
  }

  // Energy Mode
  switch (energy_mode) {
    case EnergyMode::CE:
      h5.createAttribute<std::string>("energy-mode", "continuous-energy");
      break;
    case EnergyMode::MG:
      h5.createAttribute<std::string>("energy-mode", "multi-group");
      break;
  }

  // Tracking Mode
  switch (tracking) {
    case TrackingMode::SURFACE_TRACKING:
      h5.createAttribute<std::string>("transport", "surface-tracking");
      break;
    case TrackingMode::DELTA_TRACKING:
      h5.createAttribute<std::string>("transport", "delta-tracking");
      break;
    case TrackingMode::IMPLICIT_LEAKAGE_DELTA_TRACKING:
      h5.createAttribute<std::string>("transport",
                                      "implicit-leakage-delta-tracking");
      break;
    case TrackingMode::CARTER_TRACKING:
      h5.createAttribute<std::string>("transport", "carter-tracking");
      break;
  }

  // MG CT specific
  if (energy_mode == EnergyMode::MG) {
    h5.createAttribute("ngroups", ngroups);

    h5.createAttribute("energy-bounds", energy_bounds);

    if (tracking == TrackingMode::CARTER_TRACKING) {
      h5.createAttribute("sampling-xs-ratio", sample_xs_ratio);
    }
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
  h5.createAttribute("nparticles", nparticles);

  h5.createAttribute("ngenerations", ngenerations);

  h5.createAttribute("nignored", nignored);

  h5.createAttribute("max-run-time", max_time);

  h5.createAttribute("rng-stride-warnings", rng_stride_warnings);

  h5.createAttribute("seed", rng_seed);

  h5.createAttribute("stride", rng_stride);

  h5.createAttribute("wgt-cutoff", wgt_cutoff);

  h5.createAttribute("wgt-survival", wgt_survival);

  h5.createAttribute("wgt-split", wgt_split);

  h5.createAttribute<bool>("cancellation", regional_cancellation);

  h5.createAttribute<bool>("noise-cancellation", regional_cancellation_noise);

  h5.createAttribute("cancel-noise-gens", n_cancel_noise_gens);

  h5.createAttribute<bool>("pair-distance-sqrt", pair_distance_sqrd);

  h5.createAttribute<bool>("families", families);

  h5.createAttribute<bool>("empty-entropy-bins", empty_entropy_bins);
}
}  // namespace settings
