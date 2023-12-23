/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
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
#ifndef GUI_PLOTTER_H
#define GUI_PLOTTER_H
#ifdef ABEILLE_GUI_PLOT

#include <geometry/cell.hpp>
#include <utils/direction.hpp>
#include <utils/position.hpp>

#include <ImApp/imapp.hpp>

#include <map>
#include <mutex>

namespace plotter {

class GuiPlotter : public ImApp::Layer {
 public:
  GuiPlotter();

  void render() override final;

 private:
  void render_viewport();
  void render_controls();

  ImApp::Pixel get_color(::Cell* cell);
  ImApp::Pixel get_random_color();

  Direction get_tracking_direction() const;
  Position get_start_position(uint64_t i) const;
  Direction get_comp_tracking_direction() const;
  Position get_comp_start_position(uint64_t j) const;
  Position get_pixel_position(uint64_t i, uint64_t j) const;

  void render_image();

  enum Basis : int { XY = 0, XZ = 1, YZ = 2 };
  enum ColorBy : int { Cell = 0, Material = 1 };

  // Maps for colors from plotter.hpp
  std::map<uint32_t, ImApp::Pixel> cell_id_to_color;
  std::map<uint32_t, ImApp::Pixel> material_id_to_color;
  ImApp::Image image;
  std::mutex create_color_mutex;
  int adjust_w_or_h;
  double height, width;   // Width of image in physical space ([cm]).
  double ox, oy, oz;      // Plot origin
  double mx, my, mz;      // Mouse position
  ::Cell* mcell;          // Mouse cell
  uint32_t mcell_inst;    // Mouse cell instance
  ::Material* mmaterial;  // Mouse material
  double dist_per_pixel;
  ImApp::Pixel background;
  ColorBy colorby;
  Basis basis;
  bool must_rerender;
  bool outline_boundaries;
};

}  // namespace plotter

#endif  // ABEILLE_GUI_PLOT
#endif  // GUI_PLOTTER_H
