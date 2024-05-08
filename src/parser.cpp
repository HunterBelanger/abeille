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
#include <geometry/cell_universe.hpp>
#include <geometry/geometry.hpp>
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/all_surfaces.hpp>
#include <materials/nuclide.hpp>
#include <simulation/entropy.hpp>
#include <simulation/simulation.hpp>
#include <tallies/mesh_tally.hpp>
#include <tallies/tallies.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/nd_directory.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <memory>
#include <sstream>
#include <string>

//===========================================================================
// Initialize Maps from id to index
std::map<uint32_t, size_t> surface_id_to_indx;
std::map<uint32_t, size_t> cell_id_to_indx;
std::map<uint32_t, size_t> universe_id_to_indx;

//===========================================================================
// Objects to build Simulation
std::shared_ptr<Simulation> simulation = nullptr;
std::string xspath;

void parse_input_file(std::string fname) {
  // Start parsing timer
  Timer parsing_timer;
  parsing_timer.start();

  Output::instance().write(" Reading input file.\n");

  // Open input file
  YAML::Node input = YAML::LoadFile(fname);

  if (mpi::rank == 0) {
    // Copy input file to the output file. Do this in its own scope so that
    // the potentially large amount of memory is freed before we make things.
    std::ifstream input_file_stream(fname);
    std::stringstream input_file_buffer;
    input_file_buffer << input_file_stream.rdbuf();
    std::vector<std::string> input_file_str_vec = {input_file_buffer.str()};

    auto h5 = Output::instance().h5();
    auto ds = h5.createDataSet<std::string>(
        "input", H5::DataSpace::From(input_file_str_vec));
    ds.write(input_file_str_vec);
  }

  make_settings(input);

  Output::instance().write(" Constructing simulation.\n");
  make_materials(input);

  make_geometry(input);

  make_tallies(input);

  simulation = make_simulation(input);

  // End parsing timer
  parsing_timer.stop();
  Output::instance().write(
      " Time to Parse Input : " + std::to_string(parsing_timer.elapsed_time()) +
      " seconds.\n");

  if (mpi::rank == 0) {
    Output::instance().h5().createAttribute("parse-time",
                                            parsing_timer.elapsed_time());
  }
}

void make_materials(const YAML::Node& input, bool plotting_mode) {
  // Parse materials
  if (input["materials"] && input["materials"].IsSequence()) {
    // Go through all materials
    for (size_t s = 0; s < input["materials"].size(); s++) {
      make_material(input["materials"][s], plotting_mode);
    }

  } else {
    // If materials doesn't exist and isn't a sequence, kill program
    fatal_error("No materials are provided in input file.");
  }

  // Once all materials have been read in, we need to go find the max
  // and min energy grid values. To do this, we loop through all
  // nuclides in the problem.
  if (plotting_mode) return;

  for (const auto& mat : materials) {
    double emax = mat.second->max_energy();
    double emin = mat.second->min_energy();

    if (emin > settings::min_energy) {
      settings::min_energy = emin;
    }

    if (emax < settings::max_energy) {
      settings::max_energy = emax;
    }
  }

  std::stringstream mssg;
  mssg << " Min energy = " << std::setprecision(3) << std::scientific
       << settings::min_energy << " MeV.\n";
  mssg << " Max energy = " << settings::max_energy << " MeV.\n";

  Output::instance().write(mssg.str());
}

void make_geometry(const YAML::Node& input) {
  // Parse Surfaces
  if (input["surfaces"] && input["surfaces"].IsSequence()) {
    // Go through all surfaces
    for (size_t s = 0; s < input["surfaces"].size(); s++) {
      auto surf = make_surface(input["surfaces"][s]);

      // Add surface ID to map of surface indicies
      if (surface_id_to_indx.find(surf->id()) != surface_id_to_indx.end()) {
        // ID already exists
        std::stringstream mssg;
        mssg << "The surface ID " << surf->id() << " appears multiple times.";
        fatal_error(mssg.str());
      } else {
        surface_id_to_indx[surf->id()] = geometry::surfaces.size();
        geometry::surfaces.push_back(surf);
      }
    }
  } else {
    // If surfaces doesn't exist and isn't a sequence, kill program
    fatal_error("No surfaces are provided in input file.");
  }

  // Parse Cells
  if (input["cells"] && input["cells"].IsSequence()) {
    // Go through all cells
    for (size_t c = 0; c < input["cells"].size(); c++) {
      make_cell(input["cells"][c], input);
    }

  } else {
    // If cells doesn't exist and isn't a sequence, kill program
    fatal_error("No cells are provided in input file.");
  }

  // Parse all universes
  if (input["universes"] && input["universes"].IsSequence()) {
    // Go through all cells
    for (size_t c = 0; c < input["universes"].size(); c++) {
      make_universe(input["universes"][c], input);
    }

  } else {
    // If cells doesn't exist and isn't a sequence, kill program
    fatal_error("No universes are provided in input file.");
  }

  // Make sure no universe has a cyclic dependence on itself !
  for (const auto& uni : geometry::universes) {
    if (uni->contains_universe(uni->id())) {
      std::stringstream mssg;
      mssg << "Cyclical universe reference found for universe " << uni->id()
           << ".";
      fatal_error(mssg.str());
    }
  }

  // Now that all surfaces, cells, and universes have been created, we can go
  // through and create all of the offset maps for determining the unique
  // instance of each material cell.
  for (auto& uni : geometry::universes) {
    uni->make_offset_map();
  }

  // Parse root universe
  if (input["root-universe"] && input["root-universe"].IsScalar()) {
    uint32_t root_id = input["root-universe"].as<uint32_t>();

    // Make sure it can be found
    if (universe_id_to_indx.find(root_id) == universe_id_to_indx.end()) {
      fatal_error("Root-Universe id was not found.");
    }

    geometry::root_universe = geometry::universes[universe_id_to_indx[root_id]];
  } else {
    // If doesn't exist and isn't a scalar, kill program
    fatal_error("No root-universe is provided in input file.");
  }
}

void make_universe(const YAML::Node& uni_node, const YAML::Node& input) {
  uint32_t id;
  if (uni_node["id"] && uni_node["id"].IsScalar()) {
    id = uni_node["id"].as<uint32_t>();

    // Make sure universe can't be found
    if (universe_id_to_indx.find(id) == universe_id_to_indx.end()) {
      if (uni_node["cells"] && uni_node["cells"].IsSequence()) {
        make_cell_universe(uni_node);
      } else if (uni_node["pitch"] && uni_node["pitch"].IsSequence()) {
        make_lattice(uni_node, input);
      } else {
        fatal_error("Invalid universe definition.");
      }
    }
  } else {
    fatal_error("Universe must have a valid id.");
  }
}

void find_universe(const YAML::Node& input, uint32_t id) {
  // Iterate through universes
  if (universe_id_to_indx.find(id) == universe_id_to_indx.end()) {
    bool found = false;
    for (size_t u = 0; u < input["universes"].size(); u++) {
      if (input["universes"][u]["id"] &&
          input["universes"][u]["id"].IsScalar()) {
        uint32_t u_id = input["universes"][u]["id"].as<uint32_t>();
        if (u_id == id) {
          make_universe(input["universes"][u], input);
          found = true;
          break;
        }
      } else {
        fatal_error("Universes must have a valid id.");
      }
    }

    if (found == false) {
      std::stringstream mssg;
      mssg << "Could not find universe with id " << id << ".";
      fatal_error(mssg.str());
    }
  }
}

void make_settings(const YAML::Node& input) {
  //----------------------------------------------------------------------------
  // Get the simulation mode, just for parsing input materials and making sure
  // MG data has everything needed
  if (!input["simulation"]) {
    fatal_error("No simulation entry provided in input.");
  } else if (input["simulation"] && input["simulation"].IsMap() == false) {
    fatal_error("Invalid simulation entry provided in input.");
  }
  const YAML::Node& sim = input["simulation"];

  // First, get the string identifying the simulation mode
  if (!sim["mode"]) {
    fatal_error("No mode entry in simulation definition.");
  } else if (sim["mode"].IsScalar() == false) {
    fatal_error("Invalid mode entry in simulation.");
  }
  std::string mode = sim["mode"].as<std::string>();

  if (mode == "k-eigenvalue") {
    settings::sim_mode = settings::SimMode::KEFF;
  } else if (mode == "fixed-source") {
    settings::sim_mode = settings::SimMode::FIXED_SOURCE;
  } else if (mode == "modified-fixed-source") {
    settings::sim_mode = settings::SimMode::FIXED_SOURCE;
  } else if (mode == "noise") {
    settings::sim_mode = settings::SimMode::NOISE;
  } else {
    fatal_error("Unknown simulation mode " + mode + ".");
  }

  //----------------------------------------------------------------------------
  if (input["settings"] && input["settings"].IsMap()) {
    const auto& settnode = input["settings"];

    // Get energy mode type
    if (settnode["energy-mode"] && settnode["energy-mode"].IsScalar()) {
      if (settnode["energy-mode"].as<std::string>() == "multi-group") {
        settings::energy_mode = settings::EnergyMode::MG;
        Output::instance().write(" Running in Multi-Group mode.\n");
      } else if (settnode["energy-mode"].as<std::string>() ==
                 "continuous-energy") {
        settings::energy_mode = settings::EnergyMode::CE;
        Output::instance().write(" Running in Continuous-Energy mode.\n");
      } else {
        std::string mssg = "Invalid energy mode ";
        mssg += settnode["energy-mode"].as<std::string>() + ".";
        fatal_error(mssg);
      }
    } else {
      // CE by default
      settings::energy_mode = settings::EnergyMode::CE;
      Output::instance().write(" Running in Continuous-Energy mode.\n");
    }

    // If we are in CE mode, we need to get the path to the nuclear data
    // directory file, which is either in the settings, or in the
    // environment variable. The settings file takes precedence. We also
    // need to read the temperature interpolation method too.
    if (settings::energy_mode == settings::EnergyMode::CE) {
      // Get temperature interpolation method
      if (settnode["temperature-interpolation"] &&
          settnode["temperature-interpolation"].IsScalar()) {
        std::string temp_interp =
            settnode["temperature-interpolation"].as<std::string>();

        if (temp_interp == "exact") {
          settings::temp_interpolation = TempInterpolation::Exact;
        } else if (temp_interp == "nearest") {
          settings::temp_interpolation = TempInterpolation::Nearest;
        } else if (temp_interp == "linear") {
          settings::temp_interpolation = TempInterpolation::Linear;
        } else {
          std::stringstream mssg;
          mssg << "Unknown \"temperature-interpolation\" entry in settings: "
               << temp_interp << ".";
          fatal_error(mssg.str());
        }
      } else if (settnode["temperature-interpolation"]) {
        fatal_error("Invalid \"temperature-interpolation\" entry in settings.");
      }

      // Write temperature interpolation method
      switch (settings::temp_interpolation) {
        case TempInterpolation::Exact:
          Output::instance().write(" Temperature interpolation: Exact\n");
          break;

        case TempInterpolation::Nearest:
          Output::instance().write(" Temperature interpolation: Nearest\n");
          break;

        case TempInterpolation::Linear:
          Output::instance().write(" Temperature interpolation: Linear\n");
          break;
      }

      // Get file name for nd_directory
      if (settnode["nuclear-data"] && settnode["nuclear-data"].IsScalar()) {
        // File provided in settings file
        settings::nd_directory_fname =
            settnode["nuclear-data"].as<std::string>();
      } else if (settnode["nuclear-data"]) {
        // A file entry is present, but is of invalid type.
        fatal_error("Invalid \"nulcear-data\" entry in settings.");
      } else {
        // Get the file name from the environment variable
        if (!std::getenv("ABEILLE_ND_DIRECTORY")) {
          fatal_error(
              "No \"nuclear-data\" entry in settings, and the environment "
              "variable ABEILLE_ND_DIRECTORY is undefined.");
        }
        settings::nd_directory_fname = std::getenv("ABEILLE_ND_DIRECTORY");
      }

      Output::instance().write(
          " Nuclear Data Directory: " + settings::nd_directory_fname + "\n");

      // Now we construct the global NDDirectory which lives in the settings
      settings::nd_directory = std::make_unique<NDDirectory>(
          settings::nd_directory_fname, settings::temp_interpolation);

      // Check if use-dbrc is specified
      if (settnode["use-dbrc"] && settnode["use-dbrc"].IsScalar()) {
        settings::use_dbrc = settnode["use-dbrc"].as<bool>();
      } else if (settnode["use-dbrc"]) {
        fatal_error("Invalid \"use-dbrc\" entry in settings.");
      }
      settings::nd_directory->set_use_dbrc(settings::use_dbrc);
      if (settings::use_dbrc == false) {
        Output::instance().write(" DBRC disabled in settings.\n");
      }

      // If using DBRC, check for the list of neutrons provided by the user
      if (settings::use_dbrc) {
        if (settnode["dbrc-nuclides"] &&
            settnode["dbrc-nuclides"].IsSequence()) {
          settings::dbrc_nuclides =
              settnode["dbrc-nuclides"].as<std::vector<std::string>>();
          settings::nd_directory->set_dbrc_nuclides(settings::dbrc_nuclides);
        } else if (settnode["dbrc-nuclides"]) {
          fatal_error("Invalid \"dbrc-nuclides\" entry in settings.");
        }
      }

      // Get use of URR PTables option
      if (settnode["use-urr-ptables"] &&
          settnode["use-urr-ptables"].IsScalar()) {
        settings::use_urr_ptables = settnode["use-urr-ptables"].as<bool>();
      } else if (settnode["use-urr-ptables"]) {
        fatal_error("Invalid \"use-urr-ptables\" entry in settings.");
      }
      if (settings::use_urr_ptables) {
        Output::instance().write(" Using URR PTables.\n");
      } else {
        Output::instance().write(" Not using URR PTables.\n");
      }
    }

    // If we are multi-group, get number of groups
    if (settings::energy_mode == settings::EnergyMode::MG) {
      if (!settnode["ngroups"] || !settnode["ngroups"].IsScalar()) {
        fatal_error("Number of groups for mutli-group mode not provided.");
      }
      int ngrps = settnode["ngroups"].as<int>();

      if (ngrps <= 0) {
        fatal_error("Number of groups may not be negative.");
      }

      settings::ngroups = static_cast<uint32_t>(ngrps);

      // Must also get the energy-bounds for the groups
      if (!settnode["energy-bounds"] ||
          !settnode["energy-bounds"].IsSequence()) {
        fatal_error(
            "No energy-bounds entry found in settings for multi-group mode.");
      }

      if (settnode["energy-bounds"].size() != settings::ngroups + 1) {
        fatal_error(
            "The number of energy-bounds must be equal to ngroups + 1.");
      }
      settings::energy_bounds =
          settnode["energy-bounds"].as<std::vector<double>>();

      // Check bounds
      if (!std::is_sorted(settings::energy_bounds.begin(),
                          settings::energy_bounds.end())) {
        fatal_error("The energy-bounds for multi-group mode are not sorted.");
      }

      // Make sure initial energy bound isn't exactly zero
      if (settings::energy_bounds.front() == 0.) {
        settings::energy_bounds.front() = 1.E-11;
      }
    }

    // Get roulette settings
    if (settnode["wgt-cutoff"] && settnode["wgt-cutoff"].IsScalar()) {
      settings::wgt_cutoff = settnode["wgt-cutoff"].as<double>();

      if (settings::wgt_cutoff < 0.) {
        fatal_error("Roulette cutoff (\"wgt-cutoff\") must be >= 0.");
      }
    } else if (settnode["wgt-cutoff"]) {
      fatal_error("Invalid \"wgt-cutoff\" entry in settings.");
    }

    if (settnode["wgt-survival"] && settnode["wgt-survival"].IsScalar()) {
      settings::wgt_survival = settnode["wgt-survival"].as<double>();

      if (settings::wgt_survival <= 0.) {
        fatal_error("Roulette survival weight (\"wgt-survival\") must be > 0.");
      }

      if (settings::wgt_survival <= settings::wgt_cutoff) {
        fatal_error(
            "Roulette survival weight (\"wgt-survival\") must be > cutoff "
            "weight (\"wgt-cutoff\").");
      }
    } else if (settnode["wgt-survival"]) {
      fatal_error("Invalid \"wgt-survival\" entry in settings.");
    }

    if (settnode["wgt-split"] && settnode["wgt-split"].IsScalar()) {
      settings::wgt_split = settnode["wgt-split"].as<double>();

      if (settings::wgt_split <= 0.) {
        fatal_error("Splitting weight (\"wgt-split\") must be > 0.");
      }

      if (settings::wgt_split <= settings::wgt_survival) {
        fatal_error(
            "Splitting weight (\"wgt-split\") must be > roulette survival "
            "weight (\"wgt-survival\").");
      }
    } else if (settnode["wgt-split"]) {
      fatal_error("Invalid \"wgt-split\" entry in settings.");
    }

    // Get max run time
    if (settnode["max-run-time"] && settnode["max-run-time"].IsScalar()) {
      // Get runtime in minutes
      double max_time_in_mins = settnode["max-run-time"].as<double>();
      Output::instance().write(
          " Max Run Time: " + std::to_string(max_time_in_mins) + " mins.\n");

      // Change minutes to seconds
      settings::max_time = max_time_in_mins * 60.;
    } else if (settnode["max-run-time"]) {
      fatal_error("Invalid max-run-time entry in settings.");
    }

    // See if the user wants to see RNG stride warnings
    if (settnode["rng-stride-warnings"] &&
        settnode["rng-stride-warnings"].IsScalar()) {
      settings::rng_stride_warnings =
          settnode["rng-stride-warnings"].as<bool>();
    } else if (settnode["rng-stride-warnings"]) {
      fatal_error(
          "Invalid \"rng-stride-warnings\" entry provided in settings.");
    }

    // Get seed for rng
    if (settnode["seed"] && settnode["seed"].IsScalar()) {
      settings::rng_seed = settnode["seed"].as<uint64_t>();
      Output::instance().write(
          " RNG seed = " + std::to_string(settings::rng_seed) + " ...\n");
    }

    // Get stride for rng
    if (settnode["stride"] && settnode["stride"].IsScalar()) {
      settings::rng_stride = settnode["stride"].as<uint64_t>();
      Output::instance().write(
          " RNG stride = " + std::to_string(settings::rng_seed) + " ...\n");
    }
  } else {
    // Use default settings for everything.

    // By default, we are in CE mode, so we can only get the ND directory
    // from the environment variable
    if (!std::getenv("ABEILLE_ND_DIRECTORY")) {
      fatal_error(
          "No \"nuclear-data\" entry in settings, and the environment "
          "variable ABEILLE_ND_DIRECTORY is undefined.");
    }
    settings::nd_directory_fname = std::getenv("ABEILLE_ND_DIRECTORY");

    Output::instance().write(
        " Nuclear Data Directory: " + settings::nd_directory_fname + "\n");

    // Now we construct the global NDDirectory which lives in the settings
    settings::nd_directory = std::make_unique<NDDirectory>(
        settings::nd_directory_fname, settings::temp_interpolation);
  }

  settings::write_settings_to_output();
}

void make_tallies(const YAML::Node& input) {
  // Make base tallies object which is required
  auto& tallies = Tallies::instance();

  if (!input["tallies"]) {
    return;
  } else if (!input["tallies"].IsSequence()) {
    fatal_error("Tallies entry must be provided as a sequence.");
  }

  // First, read all of the tally filters
  make_tally_filters(tallies, input);

  // Add all spatial mesh tallies to the tallies instance
  for (size_t t = 0; t < input["tallies"].size(); t++) {
    add_mesh_tally(tallies, input["tallies"][t]);
  }
}
