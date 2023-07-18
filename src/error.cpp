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
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>

void error(std::string mssg, std::source_location loc) {
  std::string message = "\n ERROR: " + mssg + "\n";
  message += " Location: " + std::string(loc.file_name()) + ":" +
             std::to_string(loc.line()) + "\n";
  Output::instance().write_error(message);
}

void fatal_error(std::string mssg, std::source_location loc) {
  std::string message = "\n FATAL ERROR: " + mssg + "\n";
  message += " Location: " + std::string(loc.file_name()) + ":" +
             std::to_string(loc.line()) + "\n";
  Output::instance().write_error(message);

  // Exit
  mpi::abort_mpi();
  std::exit(1);
}

void warning(std::string mssg, std::source_location loc) {
  std::string message = "\n WARNING: " + mssg + "\n";
  message += " Location: " + std::string(loc.file_name()) + ":" +
             std::to_string(loc.line()) + "\n";
  Output::instance().write_error(message);
}
