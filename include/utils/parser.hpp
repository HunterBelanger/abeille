/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#ifndef PARSER_H
#define PARSER_H

#include <yaml-cpp/yaml.h>

#include <map>
#include <memory>
#include <simulation/cancelator.hpp>
#include <simulation/noise_maker.hpp>
#include <simulation/noise_source.hpp>
#include <simulation/simulation.hpp>
#include <simulation/source.hpp>
#include <simulation/tallies.hpp>
#include <simulation/transporter.hpp>
#include <string>
#include <utils/settings.hpp>

//===========================================================================
// Maps to go from id to index
extern std::map<uint32_t, size_t> surface_id_to_indx;
extern std::map<uint32_t, size_t> cell_id_to_indx;
extern std::map<uint32_t, size_t> universe_id_to_indx;
extern std::map<uint32_t, size_t> lattice_id_to_indx;

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
void make_materials(YAML::Node input, bool plotting_mode = false);

// Only reads and builds geometry (used for plotting).
void make_geometry(YAML::Node input);

// Makes any type of surface from a surface yaml node
void make_surface(YAML::Node surface_node);

// Makes a universe based on wether cells or lattice
void make_universe(YAML::Node uni_node, YAML::Node input);

// Locates and then builds unknown universe
void find_universe(YAML::Node input, uint32_t id);

// Reads settings populates pointer
void make_settings(YAML::Node input);

// Reads tallies and populates pointer
void make_tallies(YAML::Node input);

// Reads into to make transporter
void make_transporter();

// Reads regional cancellation bins
void make_cancellation_bins(YAML::Node input);

// Reads sources and populates sources vector
void make_sources(YAML::Node input);

// Reads noise sources and populates noise_sources vector
void make_noise_sources(YAML::Node input);

// Construct simulation pointer
void make_simulation();

// Get the entropy mesh if given
void make_entropy_mesh(YAML::Node entropy);

#endif
