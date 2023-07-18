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
#ifndef PLOTTER_H
#define PLOTTER_H

#include <plotting/pixel.hpp>
#include <plotting/slice_plot.hpp>

#include <yaml-cpp/yaml.h>

#include <map>

namespace plotter {
// Maps for colors
extern std::map<uint32_t, Pixel> cell_id_to_color;
extern std::map<uint32_t, Pixel> material_id_to_color;

// Function to begin plotting system
void plotter(const std::string& input_fname);

#ifdef ABEILLE_GUI_PLOT
// Function to start the GUI plotter
void gui(const std::string& input_fname);
#endif

// Function to begin system to generate a slice plot
// of a selected geometry file
void slice_plotter(YAML::Node plot_node);

}  // namespace plotter

#endif
