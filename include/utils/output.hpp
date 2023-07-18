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
#ifndef OUTPUT_H
#define OUTPUT_H

#include <utils/header.hpp>
#include <utils/mpi.hpp>

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <ctime>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <source_location>
#include <string>
#include <vector>

//============================================================================
// Output Class (Singleton)
class Output {
 public:
  // No copying or moving Output
  Output(const Output&) = delete;
  Output(Output&&) = delete;
  Output& operator=(const Output&) = delete;
  Output& operator=(Output&&) = delete;
  ~Output() = default;

  static void set_output_filename(const std::string& fname);

  static Output& instance();

  template <class T>
  void write(const T& message) {
    if (mpi::rank == 0) {
      write_mutex.lock();
      std::cout << message << std::flush;
      write_mutex.unlock();
    }
  }

  void save_warning(const std::string& message);

  void write_saved_warnings();

  void write_error(const std::string& message);

  // Only master (i.e. mpi::rank == 0) can call this method without error !
  H5::File& h5(
      const std::source_location& loc = std::source_location::current());

 private:
  Output();
  static std::string output_filename;

  std::optional<H5::File> output;
  std::vector<std::string> warnings_for_latter;
  std::mutex write_mutex;
  std::mutex save_mutex;
};  // Output

//============================================================================
// Non-Member Functions for Output
void print_header();

std::string current_date_time();

#endif
