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
#ifndef PARSER_H
#define PARSER_H

#include <simulation/cancelator.hpp>
#include <simulation/noise_maker.hpp>
#include <simulation/noise_source.hpp>
#include <simulation/simulation.hpp>
#include <simulation/source.hpp>
#include <simulation/tallies.hpp>
#include <simulation/transporter.hpp>
#include <utils/settings.hpp>

#include <yaml-cpp/yaml.h>

#include <map>
#include <memory>
#include <string>

//===========================================================================
// Maps to go from id to index
extern std::map<uint32_t, size_t> surface_id_to_indx;
extern std::map<uint32_t, size_t> cell_id_to_indx;
extern std::map<uint32_t, size_t> universe_id_to_indx;

//===========================================================================
// Object to build Simulation
extern std::vector<std::shared_ptr<Source>> sources;
extern NoiseMaker noise_maker;
extern std::shared_ptr<Tallies> tallies;
extern std::shared_ptr<Transporter> transporter;
extern std::shared_ptr<Simulation> simulation;
extern std::shared_ptr<Cancelator> cancelator;
extern std::string xspath;

// Main function to parse input file
void parse_input_file(std::string fname);

// Reads all materials
void make_materials(const YAML::Node& input, bool plotting_mode = false);

// Only reads and builds geometry (used for plotting).
void make_geometry(const YAML::Node& input);

// Makes any type of surface from a surface yaml node
void make_surface(const YAML::Node& surface_node);

// Makes a universe based on wether cells or lattice
void make_universe(const YAML::Node& uni_node, const YAML::Node& input);

// Locates and then builds unknown universe
void find_universe(const YAML::Node& input, uint32_t id);

// Reads settings populates pointer
void make_settings(const YAML::Node& input);

// Reads tallies and populates pointer
void make_tallies(const YAML::Node& input);

// Reads into to make transporter
void make_transporter();

// Reads regional cancellation bins
void make_cancellation_bins(const YAML::Node& input);

// Reads sources and populates sources vector
void make_sources(const YAML::Node& input);

// Reads noise sources and populates noise_sources vector
void make_noise_sources(const YAML::Node& input);

// Construct simulation pointer
void make_simulation();

// Get the entropy mesh if given
void make_entropy_mesh(const YAML::Node& entropy);

#endif
