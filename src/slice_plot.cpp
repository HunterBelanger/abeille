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
#include <geometry/cell.hpp>
#include <geometry/geometry.hpp>
#include <plotting/slice_plot.hpp>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <fstream>
#include <memory>

namespace plotter {

SlicePlot::SlicePlot(std::string fname, uint64_t pwidth, uint64_t pheight,
                     double width, double height, Position origin, Basis basis,
                     ColorBy colorby, Pixel background)
    : file_name_{fname},
      plot_width_{pwidth},
      plot_height_{pheight},
      basis_{basis},
      colorby_{colorby},
      origin_{origin},
      background_{background},
      image_matrix(),
      width_{width},
      height_{height},
      create_color_mutex{} {}

Pixel SlicePlot::get_color(::Cell* cell) {
  Pixel pixel_color = background_;
  // If pointer isn't nullpntr, get pixel
  if (cell != nullptr) {
    if (colorby_ == ColorBy::Cell) {
      // Do same check twice with mutex to make thread safe
      if (cell_id_to_color.find(cell->id()) == cell_id_to_color.end()) {
        // Check if cell id is in id_to_pixel
        create_color_mutex.lock();
        if (cell_id_to_color.find(cell->id()) == cell_id_to_color.end()) {
          // Get new random color for id
          cell_id_to_color[cell->id()] = get_random_color();
        }
        create_color_mutex.unlock();
      }
      pixel_color = cell_id_to_color[cell->id()];
    } else {
      // Color by material
      const ::Material* material = cell->material();
      // Do same check twice with mutex to make thread safe
      if (material_id_to_color.find(material->id()) ==
          material_id_to_color.end()) {
        // Check if cell id is in id_to_pixel
        create_color_mutex.lock();
        if (material_id_to_color.find(material->id()) ==
            material_id_to_color.end()) {
          // Get new random color for id
          material_id_to_color[material->id()] = get_random_color();
        }
        create_color_mutex.unlock();
      }
      pixel_color = material_id_to_color[material->id()];
    }
  }
  return pixel_color;
}

Direction SlicePlot::get_tracking_direction() const {
  switch (basis_) {
    case Basis::XY:
      return Direction(0., 1., 0.);
      break;
    case Basis::YZ:
      return Direction(0., 0., 1.);
      break;
    case Basis::XZ:
      return Direction(0., 0., 1.);
      break;
  }

  // NEVER GETS HERE
  return Direction(1., 0., 0.);
}

Position SlicePlot::get_start_position(uint64_t i) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= plot_height_) {
    fatal_error("Trying to deffine invalid pixel for plot.");
  }

  if (basis_ == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double y_low = origin_.y() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double y = y_low;

    // Return coordinate
    return {x, y, origin_.z()};
  } else if (basis_ == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    double dy = height_ / static_cast<double>(plot_height_);

    // Get upper corner off plot
    double y_low = origin_.y() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double y = (static_cast<double>(i) + 0.5) * dy + y_low;
    double z = z_low;

    // Return coordinate
    return {origin_.x(), y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double z = z_low;

    // Return coordinate
    return {x, origin_.y(), z};
  }
}

void SlicePlot::generate_plot() {
  // Calculate number of elements
  uint64_t n_pixels = plot_width_ * plot_height_;

  // Shape image matrix accordingly
  image_matrix.resize(n_pixels);

  // Get the tracking direction
  const Direction u = this->get_tracking_direction();
  const double dist_per_pixel = width_ / static_cast<double>(plot_width_);

// Go through each pixel
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < plot_height_; i++) {
    Position strt = get_start_position(i);

    // Initalize the tracker
    Tracker trkr(strt, u);
    ::Cell* cell = trkr.cell();
    Pixel pixel_color = this->get_color(cell);

    uint64_t j = 0;
    while (j < plot_width_) {
      // Get the boundary
      auto bound = trkr.get_nearest_boundary();

      // Get the number of pixels till the boundary
      const double pixels_to_bound = bound.distance / dist_per_pixel;
      uint64_t npixels = static_cast<uint64_t>(std::round(pixels_to_bound));
      if (npixels > (plot_width_ - j) || bound.distance == INF) {
        npixels = plot_width_ - j;
      }
      const double npixels_dist = static_cast<double>(npixels) * dist_per_pixel;

      // Set all pixels
      for (uint64_t p = 0; p < npixels; p++) {
        uint64_t indx = i * plot_width_ + j;
        image_matrix.at(indx) = pixel_color;
        j++;
      }

      // Cross boundary, update cell, and get new pixel
      if (pixels_to_bound - static_cast<double>(npixels) < 0.5) {
        trkr.move(npixels_dist + dist_per_pixel);
      } else {
        trkr.move(npixels_dist);
      }
      trkr.restart_get_current();
      cell = trkr.cell();
      pixel_color = this->get_color(cell);
      if (pixels_to_bound - static_cast<double>(npixels) < 0.5) {
        uint64_t indx = i * plot_width_ + j;
        image_matrix.at(indx) = pixel_color;
        j++;
      }
    }  // while j < plot_width_
  }  // For i which is parallel
}

void SlicePlot::write() const {
  if (image_matrix.size() > 0) {
    std::ofstream file(file_name_ + ".ppm");
    file << "P6\n";

    // Get width and height as string
    std::string width = std::to_string(plot_width_);
    std::string height = std::to_string(plot_height_);
    file << width << ' ' << height << '\n';

    // Set max value
    file << "255\n";

    // Write all pixels from matrix. First get pointer
    const char* img_data = reinterpret_cast<const char*>(image_matrix.data());
    // Then get number of bytes. Use sizeof to be sure works on all systems
    std::streamsize nbytes =
        static_cast<std::streamsize>(image_matrix.size() * sizeof(Pixel));
    // Write instead of streaming seems to be much faster
    file.write(img_data, nbytes);

    // Close file
    file.close();
  }
}

Position SlicePlot::get_pixel_position(uint64_t i, uint64_t j) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= plot_height_ || j >= plot_width_) {
    fatal_error("Trying to deffine invalid pixel for plot.");
  }

  if (basis_ == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);
    double dy = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double y_low = origin_.y() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double y = (static_cast<double>(j) + 0.5) * dy + y_low;

    // Return coordinate
    return {x, y, origin_.z()};
  } else if (basis_ == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    double dy = height_ / static_cast<double>(plot_height_);
    double dz = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double y_low = origin_.y() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double y = (static_cast<double>(i) + 0.5) * dy + y_low;
    double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {origin_.x(), y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);
    double dz = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {x, origin_.y(), z};
  }
}

Pixel SlicePlot::get_random_color() {
  uint8_t r = static_cast<uint8_t>(255.0 * settings::rng());
  uint8_t g = static_cast<uint8_t>(255.0 * settings::rng());
  uint8_t b = static_cast<uint8_t>(255.0 * settings::rng());

  return {r, g, b};
}

}  // namespace plotter
