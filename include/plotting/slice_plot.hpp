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
#ifndef SLICE_PLOT_H
#define SLICE_PLOT_H

#include <geometry/cell.hpp>
#include <plotting/pixel.hpp>
#include <utils/direction.hpp>
#include <utils/position.hpp>

#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace plotter {
// Maps for colors from plotter.hpp
extern std::map<uint32_t, Pixel> cell_id_to_color;
extern std::map<uint32_t, Pixel> material_id_to_color;

class SlicePlot {
 public:
  enum Basis { XY, XZ, YZ };
  enum ColorBy { Cell, Material };

  SlicePlot(std::string fname, uint64_t pwidth, uint64_t pheight, double width,
            double height, Position origin, Basis basis, ColorBy colorby,
            Pixel background = Pixel(255, 255, 255));
  ~SlicePlot() = default;

  void generate_plot();
  void write() const;

 private:
  std::string file_name_;
  uint64_t plot_width_, plot_height_;  // Both are in number of pixels
  Basis basis_;
  ColorBy colorby_;
  Position origin_;
  Pixel background_;
  std::vector<Pixel> image_matrix;  // Row-Major order
  double width_, height_;
  std::mutex create_color_mutex;

  Position get_pixel_position(uint64_t i, uint64_t j) const;
  Position get_start_position(uint64_t i) const;
  Direction get_tracking_direction() const;
  Pixel get_color(::Cell* cell);
  Pixel get_random_color();

};  // SlicePlot

}  // namespace plotter

#endif
