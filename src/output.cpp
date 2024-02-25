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
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

#ifdef ABEILLE_USE_OMP
#include <omp.h>  // For print_header function
#endif

//============================================================================
// Initialization of static members of Output singleton
std::string Output::output_filename = "miel.h5";

//============================================================================
// Output Singleton Methods
Output::Output() : output(std::nullopt), warnings_for_latter() {
  // Only open the output file is we are master
  if (mpi::rank == 0 && settings::plotting_mode == false) {
    // First, we see if the output_filename already exists. If it does, we
    // delete it.
    if (std::filesystem::exists(output_filename)) {
      bool did_rm = std::filesystem::remove(output_filename);
      if (did_rm == false) {
        // Cannot use fatal_error here, because that calls Output !
        std::cerr
            << "\n FATAL ERROR: Could not delete pre-existing file/directory "
            << output_filename << ".\n";
        std::cerr << std::flush;
        std::exit(1);
      }
    }

    try {
      output = H5::File(output_filename, H5::File::Create);
    } catch (...) {
      // Cannot use fatal_error here, because that calls Output !
      std::cerr << "\n FATAL ERROR: Error creating output file "
                << output_filename << ".\n";
      std::cerr << std::flush;
      std::exit(1);
    }
  }
}

void Output::set_output_filename(const std::string& fname) {
  output_filename = fname;
}

Output& Output::instance() {
  static Output output_instance;

  return output_instance;
}

void Output::save_warning(const std::string& message) {
  save_mutex.lock();
  warnings_for_latter.push_back(message);
  save_mutex.unlock();
}

void Output::write_saved_warnings() {
  if (!warnings_for_latter.empty()) {
    std::cout << "\n";
  }

  for (const auto& wrn : warnings_for_latter) {
    std::cout << " WARNING: " << wrn << "\n";
  }

  if (!warnings_for_latter.empty()) {
    std::cout << "\n";
  }

  warnings_for_latter.clear();
}

void Output::write_error(const std::string& message) {
  write_mutex.lock();
  std::cerr << message << std::flush;
  write_mutex.unlock();
}

H5::File& Output::h5() {
  if (mpi::rank != 0) {
    std::cerr << " FATAL ERROR: Only master can write to output file.\n";
    std::exit(1);
    fatal_error("Only master can write to output file.");
  }

  return *output;
}

//============================================================================
// Non-Member Functions
void print_header() {
  Output& output = Output::instance();
  output.write(logo);
  output.write(header);

  // Version info
  std::string info = "";
  info += " Version             : " + std::string(ABEILLE_VERSION_STRING)
#if defined(DEVELOPMENT_VERSION)
          + " (Development)"
#endif
          + "\n";

  // Git Hash
  info += " Abeille Git Hash    : " ABEILLE_GIT_HASH "\n";

  // Build date
  info += " Build Date          : " + std::string(__DATE__) + " ";
  info += std::string(__TIME__) + "\n";

  // Current date and time
  info += " Date/Time           : " + current_date_time() + "\n";

#ifdef ABEILLE_USE_OMP
  info +=
      " OpenMP Threads      : " + std::to_string(omp_get_max_threads()) + "\n";
#endif

#ifdef ABEILLE_USE_MPI
  info += " MPI Ranks           : " + std::to_string(mpi::size) + "\n";
#endif

  output.write(info);
}

std::string current_date_time() {
  std::time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%b %e %Y %H:%M:%S", &tstruct);
  return buf;
}
