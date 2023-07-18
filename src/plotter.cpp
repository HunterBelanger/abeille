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
#include <plotting/plotter.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>

#ifdef ABEILLE_GUI_PLOT
#include <ImApp/imapp.hpp>
#include <plotting/gui_plotter.hpp>
#include "logo.cpp"
#endif

namespace plotter {
// Maps for plotting colors
std::map<uint32_t, Pixel> cell_id_to_color;
std::map<uint32_t, Pixel> material_id_to_color;

void plotter(const std::string& input_fname) {
  // Get output pointer
  Output::instance().write(" Starting plotting engine...\n");

  // Open the YAML node for input file
  YAML::Node input = YAML::LoadFile(input_fname);

  // If there are plots, load geometry
  if (input["plots"] && input["plots"].IsSequence()) {
    // Load materials in plotting_mode = true
    make_materials(input, true);

    // Load geometry portions of input
    make_geometry(input);

    // Go through all plots
    for (size_t p = 0; p < input["plots"].size(); p++) {
      if (input["plots"][p]["type"] && input["plots"][p]["type"].IsScalar()) {
        std::string plot_type = input["plots"][p]["type"].as<std::string>();
        if (plot_type == "slice") {
          slice_plotter(input["plots"][p]);
        } else {
          error(plot_type + " is not a valid plot type.\n");
        }
      } else {
        error(" Plot " + std::to_string(p) +
              " is missing a valid type attribute.\n");
      }
    }
  } else {
    error(" No plots are specified in the input file.\n");
  }
}

void slice_plotter(YAML::Node plot_node) {
  // Build SlicePlot from node
  // Get name
  std::string name = "";
  if (plot_node["name"] && plot_node["name"].IsScalar()) {
    name = plot_node["name"].as<std::string>();
  } else {
    fatal_error("Plot definition is missing name.");
  }

  // Get basis
  SlicePlot::Basis basis = SlicePlot::Basis::XY;
  std::string basis_str;
  if (plot_node["basis"] && plot_node["basis"].IsScalar()) {
    basis_str = plot_node["basis"].as<std::string>();
    if (basis_str == "xy")
      basis = SlicePlot::Basis::XY;
    else if (basis_str == "yz")
      basis = SlicePlot::Basis::YZ;
    else if (basis_str == "xz")
      basis = SlicePlot::Basis::XZ;
    else {
      fatal_error("Plot definition has an invalid basis");
    }
  } else {
    fatal_error("Plot definition is missing basis.");
  }

  // Get resolution
  uint64_t pwidth = 0, pheight = 0;
  if (plot_node["resolution"] && plot_node["resolution"].IsSequence()) {
    if (plot_node["resolution"].size() == 2) {
      pheight = plot_node["resolution"][0].as<uint64_t>();
      pwidth = plot_node["resolution"][1].as<uint64_t>();
    } else {
      fatal_error("Plot resolution must have two entries.");
    }
  } else {
    fatal_error("Plot must have a valid resolution entry.");
  }

  // Get dimensions
  double height = 0, width = 0;
  if (plot_node["dimensions"] && plot_node["dimensions"].IsSequence()) {
    if (plot_node["dimensions"].size() == 2) {
      height = plot_node["dimensions"][0].as<double>();
      width = plot_node["dimensions"][1].as<double>();
    } else {
      fatal_error("Plot dimension must have two entries.");
    }
  } else {
    fatal_error("Plot must have a valid dimensions entry.");
  }

  if (height <= 0. || width <= 0.) {
    fatal_error("Plot heigh and width must be greater than zero.");
  }

  // Get origin
  double x = 0, y = 0, z = 0;
  if (plot_node["origin"] && plot_node["origin"].IsSequence()) {
    if (plot_node["origin"].size() == 3) {
      x = plot_node["origin"][0].as<double>();
      y = plot_node["origin"][1].as<double>();
      z = plot_node["origin"][2].as<double>();
    } else {
      fatal_error("Plot origin must have three entries.");
    }
  } else {
    fatal_error("Plot must have a valid origin entry.");
  }
  Position origin(x, y, z);

  // Get Color scheme
  SlicePlot::ColorBy color = SlicePlot::ColorBy::Material;
  if (plot_node["color"] && plot_node["color"].IsScalar()) {
    if (plot_node["color"].as<std::string>() == "cell") {
      color = SlicePlot::ColorBy::Cell;
    } else if (plot_node["color"].as<std::string>() == "material") {
      color = SlicePlot::ColorBy::Material;
    } else {
      fatal_error("Plot must have a valid color scheme.");
    }
  } else {
    fatal_error("Plot must have a valid color entry.");
  }

  // Make plot
  SlicePlot plot(name, pwidth, pheight, width, height, origin, basis, color);

  Output::instance().write(" Generating " + name + " plot...\n");
  plot.generate_plot();

  plot.write();
}

#ifdef ABEILLE_GUI_PLOT
void gui(const std::string& input_fname) {
  // Get output pointer
  Output::instance().write(" Starting GUI plotter...\n");

  // Open the YAML node for input file
  YAML::Node input = YAML::LoadFile(input_fname);

  // Load materials in plotting_mode = true
  make_materials(input, true);

  // Load geometry portions of input
  make_geometry(input);

  // Build logo image
  ImApp::Image logo(LOGO_HEIGHT, LOGO_WIDTH);
  const ImApp::Pixel* pixels =
      reinterpret_cast<const ImApp::Pixel*>(&LOGO_DATA[0]);
  for (std::size_t i = 0; i < logo.size(); i++) {
    logo[i] = pixels[i];
  }

  try {
    ImApp::App guiplotter(1920, 1080, "Abeille Geometry Plotter");
    // Viewports don't always work well on Linux.
    // Might want to disable this at some point.
    guiplotter.enable_viewports();
    guiplotter.enable_docking();
    guiplotter.push_layer(std::make_unique<GuiPlotter>());
    guiplotter.set_icon(logo);
    guiplotter.run();
  } catch (std::exception& error) {
    Output::instance().write_error(error.what());
    std::exit(1);
  }
}
#endif

}  // namespace plotter
