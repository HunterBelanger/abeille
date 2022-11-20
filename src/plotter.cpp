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
#include <plotting/plotter.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>

namespace plotter {
// Maps for plotting colors
std::map<uint32_t, Pixel> cell_id_to_color;
std::map<uint32_t, Pixel> material_id_to_color;

void plotter(std::string input_fname) {
  // Get output pointer
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Starting plotting engine...\n");

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
          error(plot_type + " is not a valid plot type.\n", __FILE__, __LINE__);
        }
      } else {
        error(" Plot " + std::to_string(p) +
                  " is missing a valid type attribute.\n",
              __FILE__, __LINE__);
      }
    }
  } else {
    error(" No plots are specified in the input file.\n", __FILE__, __LINE__);
  }
}

void slice_plotter(YAML::Node plot_node) {
  // Build SlicePlot from node
  // Get name
  std::string name = "";
  if (plot_node["name"] && plot_node["name"].IsScalar()) {
    name = plot_node["name"].as<std::string>();
  } else {
    std::string mssg = "Plot definition is missing name.";
    fatal_error(mssg, __FILE__, __LINE__);
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
      std::string mssg = "Plot definition has an invalid basis";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg = "Plot definition is missing basis.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get resolution
  uint64_t pwidth = 0, pheight = 0;
  if (plot_node["resolution"] && plot_node["resolution"].IsSequence()) {
    if (plot_node["resolution"].size() == 2) {
      pheight = plot_node["resolution"][0].as<uint64_t>();
      pwidth = plot_node["resolution"][1].as<uint64_t>();
    } else {
      std::string mssg = "Plot resolution must have two entries.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg = "Plot must have a valid resolution entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get dimensions
  double height = 0, width = 0;
  if (plot_node["dimensions"] && plot_node["dimensions"].IsSequence()) {
    if (plot_node["dimensions"].size() == 2) {
      height = plot_node["dimensions"][0].as<double>();
      width = plot_node["dimensions"][1].as<double>();
    } else {
      std::string mssg = "Plot dimension must have two entries.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg = "Plot must have a valid dimensions entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (height <= 0. || width <= 0.) {
    std::string mssg = "Plot heigh and width must be greater than zero.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get origin
  double x = 0, y = 0, z = 0;
  if (plot_node["origin"] && plot_node["origin"].IsSequence()) {
    if (plot_node["origin"].size() == 3) {
      x = plot_node["origin"][0].as<double>();
      y = plot_node["origin"][1].as<double>();
      z = plot_node["origin"][2].as<double>();
    } else {
      std::string mssg = "Plot origin must have three entries.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg = "Plot must have a valid origin entry.";
    fatal_error(mssg, __FILE__, __LINE__);
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
      std::string mssg = "Plot must have a valid color scheme.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    std::string mssg = "Plot must have a valid color entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make plot
  SlicePlot plot(name, pwidth, pheight, width, height, origin, basis, color);

  Output::instance()->write(" Generating " + name + " plot...\n");
  plot.generate_plot();

  plot.write();
}

}  // namespace plotter
