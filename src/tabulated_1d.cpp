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
#include <utils/tabulated_1d.hpp>

#include <sstream>
#include <vector>

pndl::Tabulated1D make_tabulated_1d(const YAML::Node& node,
                                    const std::string& x_name,
                                    const std::string& y_name) {
  if (!node[x_name] || node[x_name].IsSequence() == false) {
    std::stringstream mssg;
    mssg << "Function node does not have a valid x-axis entry by name of \""
         << x_name << "\".";
    fatal_error(mssg.str());
  }

  std::vector<double> x_vals = node[x_name].as<std::vector<double>>();
  if (x_vals.size() < 2) {
    fatal_error("Must have at least 2 tabulated points for function.");
  }

  if (!node[y_name] || node[y_name].IsSequence() == false) {
    std::stringstream mssg;
    mssg << "Function node does not have a valid y-axis entry by name of \""
         << y_name << "\".";
    fatal_error(mssg.str());
  }

  std::vector<double> y_vals = node[y_name].as<std::vector<double>>();
  if (y_vals.size() != x_vals.size()) {
    fatal_error("Different number of x and y values.");
  }

  return pndl::Tabulated1D(pndl::Interpolation::LinLin, x_vals, y_vals);
}