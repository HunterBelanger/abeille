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
#ifndef SETTINGS_H
#define SETTINGS_H

#include <utils/nd_directory.hpp>
#include <utils/timer.hpp>

#include <pcg_random.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using TempInterpolation = NDDirectory::TemperatureInterpolation;

namespace settings {
enum class SimulationMode {
  K_EIGENVALUE,
  FIXED_SOURCE,
  MODIFIED_FIXED_SOURCE,
  NOISE
};
enum class TrackingMode {
  SURFACE_TRACKING,
  DELTA_TRACKING,
  IMPLICIT_LEAKAGE_DELTA_TRACKING,
  CARTER_TRACKING
};
enum class CollisionMode {
  BRANCHING,
  BRANCHLESS_ISOTOPE,
  BRANCHLESS_MATERIAL
};
enum class EnergyMode { CE, MG };

extern int nparticles;
extern int ngenerations;
extern int nignored;
extern int nskip;
extern uint32_t ngroups;
extern int n_cancel_noise_gens;

extern bool plotting_mode;

extern Timer alpha_omega_timer;
extern double max_time;

extern double min_energy;
extern double max_energy;
extern double target_at_rest_threshold;

extern SimulationMode mode;
extern TrackingMode tracking;
extern CollisionMode collisions;
extern EnergyMode energy_mode;

extern uint64_t rng_seed;
extern uint64_t rng_stride;
extern pcg32 rng;

extern double wgt_cutoff;
extern double wgt_survival;
extern double wgt_split;

extern double w_noise;
extern double eta;
extern double keff;

extern bool use_urr_ptables;

extern bool converged;

extern bool pair_distance_sqrd;
extern bool families;
extern bool empty_entropy_bins;

extern bool regional_cancellation;
extern bool regional_cancellation_noise;
extern bool inner_generations;
extern bool normalize_noise_source;
extern bool rng_stride_warnings;
extern bool load_source_file;

// PI settings
extern bool combing;
extern bool branchless_splitting;

// Energy bounds for multi-group mode
extern std::vector<double> energy_bounds;
std::size_t group(double E);
extern std::vector<double> sample_xs_ratio;
extern bool chi_matrix;
extern bool use_virtual_collisions;

extern std::string output_file_name;
extern std::string source_file_name;
extern std::string in_source_file_name;

extern std::unique_ptr<NDDirectory> nd_directory;
extern std::string nd_directory_fname;
extern TempInterpolation temp_interpolation;
extern bool use_dbrc;
extern std::vector<std::string> dbrc_nuclides;
void initialize_nd_directory();

void initialize_global_rng();

void write_settings_to_output();
}  // namespace settings

#endif  // MG_SETTINGS_H
