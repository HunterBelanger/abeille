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
#include <utils/rng.hpp>
#include <utils/timer.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using TempInterpolation = NDDirectory::TemperatureInterpolation;

namespace settings {

// We at least need a global mode entry, so that we can ensure MG data is
// complete for the simulation mode
enum class SimMode { KEFF, FIXED_SOURCE, NOISE, ALPHA };
extern SimMode sim_mode;

extern bool plotting_mode;

extern Timer alpha_omega_timer;
extern double max_time;

// Physics settings
extern double min_energy;
extern double max_energy;
extern double target_at_rest_threshold;
extern bool use_urr_ptables;

// Random Number Generator settings
extern uint64_t rng_seed;
extern uint64_t rng_stride;
extern RNG rng;
extern bool rng_stride_warnings;

// Roulette and Splitting parameters
extern double wgt_cutoff;
extern double wgt_survival;
extern double wgt_split;

enum class EnergyMode { CE, MG };
extern EnergyMode energy_mode;
// Energy bounds for multi-group mode
extern std::vector<double> energy_bounds;
extern uint32_t ngroups;
std::size_t group(double E);

extern bool chi_matrix;
extern bool use_virtual_collisions;

// Nuclear Data Directory Settings
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
