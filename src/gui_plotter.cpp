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
#ifdef ABEILLE_GUI_PLOT

#include <plotting/gui_plotter.hpp>
#include <plotting/pixel.hpp>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <filesystem>

namespace plotter {

// Maps for colors from plotter.hpp
extern std::map<uint32_t, Pixel> cell_id_to_color;
extern std::map<uint32_t, Pixel> material_id_to_color;

GuiPlotter::GuiPlotter()
    : cell_id_to_color(),
      material_id_to_color(),
      image(500, 500),
      create_color_mutex(),
      height(10.),
      width(10.),
      ox(0.),
      oy(0.),
      oz(0.),
      mx(0.),
      my(0.),
      mz(0.),
      mcell(nullptr),
      mmaterial(nullptr),
      background(),
      colorby(ColorBy::Material),
      basis(Basis::XY),
      must_rerender(true),
      outline_boundaries(true) {
  // First, fill locatl color maps with the outer color maps
  for (const auto& color : plotter::cell_id_to_color) {
    const uint32_t& id = color.first;
    const Pixel& p = color.second;
    this->cell_id_to_color[id] = ImApp::Pixel(p.R(), p.G(), p.B());
  }

  for (const auto& color : plotter::material_id_to_color) {
    const uint32_t& id = color.first;
    const Pixel& p = color.second;
    this->material_id_to_color[id] = ImApp::Pixel(p.R(), p.G(), p.B());
  }
}

void GuiPlotter::render() {
  // We put a Dockspace over the entire viewport.
  ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

  this->render_viewport();
  this->render_controls();
}

void GuiPlotter::render_viewport() {
  // Get IO instance, as we will need it for certain things
  const auto& io = ImGui::GetIO();

  // Make window
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Viewport");

  // First, get the window size. If the size doesn't match the current
  // image size, we must set must_rerender to true.
  ImVec2 size = ImGui::GetContentRegionAvail();
  std::uint32_t wwidth = static_cast<std::uint32_t>(size[0]);
  std::uint32_t wheight = static_cast<std::uint32_t>(size[1]);
  if ((wwidth != image.width()) || (wheight != image.height())) {
    image.resize(wheight, wwidth);

    dist_per_pixel = width / static_cast<double>(image.width());
    height = static_cast<double>(image.height()) * dist_per_pixel;

    must_rerender = true;
  }

  // Now we check to see if the user is dragging their mouse, to
  // change the origin of the plot
  if (ImGui::IsWindowHovered() && ImGui::IsMouseDragging(0)) {
    ImVec2 mouse_drag = io.MouseDelta;
    if (mouse_drag[0] != 0. || mouse_drag[1] != 0.) {
      switch (basis) {
        case Basis::XY:
          ox -= dist_per_pixel * mouse_drag[1];
          oy -= dist_per_pixel * mouse_drag[0];
          break;

        case Basis::XZ:
          ox -= dist_per_pixel * mouse_drag[1];
          oz -= dist_per_pixel * mouse_drag[0];
          break;

        case Basis::YZ:
          oy -= dist_per_pixel * mouse_drag[1];
          oz -= dist_per_pixel * mouse_drag[0];
          break;
      }
      must_rerender = true;
    }
  }

  // Now we check to see if the user is scrolling, which can
  // resize the image zoom, by chaning the width and height.
  if (ImGui::IsWindowHovered() && std::abs(io.MouseWheel) > 0.1) {
    // Change width
    width += 0.05 * io.MouseWheel * width;
    if (width < 1.E-6) width = 1.E-6;

    // Must recalculate height
    dist_per_pixel = width / static_cast<double>(image.width());
    height = static_cast<double>(image.height()) * dist_per_pixel;

    must_rerender = true;
  }

  // If we window must be rerendered, we do that now
  if (must_rerender) {
    this->render_image();
    must_rerender = false;

    // Now we need to send the image to the GPU.
    image.send_to_gpu();
  }

  // Get upper-left corner of image, in Window frame (pixels)
  const ImVec2 img_pos = ImGui::GetCursorPos();

  // Now we need to add the image to the window
  void* texture_id = reinterpret_cast<void*>(
      static_cast<intptr_t>(image.ogl_texture_id().value()));
  ImGui::Image(texture_id, ImVec2(static_cast<float>(image.width()),
                                  static_cast<float>(image.height())));

  // Now we get the mouse position, so the user can identify cells
  // and materials
  if (ImGui::IsWindowHovered() && !ImGui::IsMouseDragging(0)) {
    // Get the new mouse coordinates in screen space (pixels)
    const ImVec2 mouse_pos = ImGui::GetMousePos();

    // Get window position in screen space (pixels)
    const ImVec2 window_pos = ImGui::GetWindowPos();

    // Get image in screen space (pixel)
    ImVec2 img_pos_on_screen;
    img_pos_on_screen[0] = img_pos[0] + window_pos[0];
    img_pos_on_screen[1] = img_pos[1] + window_pos[1];

    // Get the position of mouse relative to image (pixel)
    ImVec2 mouse_img_pos;
    mouse_img_pos[0] = mouse_pos[0] - img_pos_on_screen[0];
    mouse_img_pos[1] = mouse_pos[1] - img_pos_on_screen[1];

    mouse_img_pos[0] -= 0.5f * static_cast<float>(image.width());
    mouse_img_pos[1] -= 0.5f * static_cast<float>(image.height());

    // Convert the image position to physical position
    switch (basis) {
      case Basis::XY:
        mx = ox + dist_per_pixel * mouse_img_pos[1];
        my = oy + dist_per_pixel * mouse_img_pos[0];
        mz = oz;
        break;

      case Basis::XZ:
        mx = ox + dist_per_pixel * mouse_img_pos[1];
        mz = oz + dist_per_pixel * mouse_img_pos[0];
        my = oy;
        break;

      case Basis::YZ:
        my = oy + dist_per_pixel * mouse_img_pos[1];
        mz = oz + dist_per_pixel * mouse_img_pos[0];
        mx = ox;
        break;
    }

    // Initialize a tracker with mouse position
    Position mp(mx, my, mz);
    Direction mu(1., 0., 0.);
    Tracker mtrkr(mp, mu);

    if (mtrkr.is_lost()) {
      mcell = nullptr;
      mcell_inst = 0;
      mmaterial = nullptr;
    } else {
      mcell = mtrkr.cell();
      mcell_inst = mtrkr.cell_instance();
      mmaterial = mtrkr.material();
    }
  }

  // Capture a right-click to bring up color changer
  if (ImGui::IsWindowHovered() && ImGui::IsMouseClicked(1) && mcell)
    ImGui::OpenPopup("Select Color");
  if (ImGui::BeginPopup("Select Color")) {
    ImApp::Pixel color = get_color(mcell);
    ImVec4 fcolor(static_cast<float>(color.r()) / 255.f,
                  static_cast<float>(color.g()) / 255.f,
                  static_cast<float>(color.b()) / 255.f,
                  static_cast<float>(color.a()) / 255.f);

    if (colorby == ColorBy::Material) {
      ImGui::Text("Material ID: %i", mmaterial->id());
      ImGui::Text("Material Name: %s", mmaterial->name().data());
    } else {
      ImGui::Text("Cell ID: %i", mcell->id());
      ImGui::Text("Cell Name: %s", mcell->name().data());
    }

    if (ImGui::ColorEdit3("", reinterpret_cast<float*>(&fcolor))) {
      color.r() = static_cast<uint8_t>(fcolor.x * 255.f);
      color.g() = static_cast<uint8_t>(fcolor.y * 255.f);
      color.b() = static_cast<uint8_t>(fcolor.z * 255.f);

      if (colorby == ColorBy::Material)
        material_id_to_color[mmaterial->id()] = color;
      else
        cell_id_to_color[mcell->id()] = color;

      must_rerender = true;
    }
    ImGui::EndPopup();
  }

  ImGui::End();
}

void GuiPlotter::render_controls() {
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Controls");

  // Origin
  ImGui::Separator();
  ImGui::Text("Plot Origin");
  if (ImGui::InputDouble("X [cm]", &ox, 0., 0.)) must_rerender = true;
  if (ImGui::InputDouble("Y [cm]", &oy, 0., 0.)) must_rerender = true;
  if (ImGui::InputDouble("Z [cm]", &oz, 0., 0.)) must_rerender = true;

  // Physical Dimensions of plot
  ImGui::Separator();
  ImGui::Text("Width/Height");
  ImGui::RadioButton("Width", &adjust_w_or_h, 0);
  ImGui::SameLine();
  ImGui::RadioButton("Height", &adjust_w_or_h, 1);
  if (adjust_w_or_h == 0) {
    if (ImGui::InputDouble("Width [cm]", &width, 0., 0.)) {
      must_rerender = true;

      if (width < 1.E-6) width = 1.E-6;

      // Must recalculate height
      dist_per_pixel = width / static_cast<double>(image.width());
      height = static_cast<double>(image.height()) * dist_per_pixel;
    }
  } else if (adjust_w_or_h == 1) {
    if (ImGui::InputDouble("Height [cm]", &height, 0., 0.)) {
      must_rerender = true;

      if (height < 1.E-6) height = 1.E-6;

      // Must recalculate width
      dist_per_pixel = height / static_cast<double>(image.height());
      width = static_cast<double>(image.width()) * dist_per_pixel;
    }
  }
  ImGui::Text("Width [cm]: %f, Height [cm]: %f", width, height);

  // Slice Basis
  ImGui::Separator();
  ImGui::Text("Slice Basis");
  if (ImGui::RadioButton("XY", reinterpret_cast<int*>(&basis), Basis::XY))
    must_rerender = true;
  ImGui::SameLine();
  if (ImGui::RadioButton("XZ", reinterpret_cast<int*>(&basis), Basis::XZ))
    must_rerender = true;
  ImGui::SameLine();
  if (ImGui::RadioButton("YZ", reinterpret_cast<int*>(&basis), Basis::YZ))
    must_rerender = true;

  // Color Method
  ImGui::Separator();
  ImGui::Text("Coloring");
  if (ImGui::RadioButton("Cell", reinterpret_cast<int*>(&colorby),
                         ColorBy::Cell))
    must_rerender = true;
  ImGui::SameLine();
  if (ImGui::RadioButton("Material", reinterpret_cast<int*>(&colorby),
                         ColorBy::Material))
    must_rerender = true;

  if (ImGui::Checkbox("Mark Boundaries", &outline_boundaries))
    must_rerender = true;

  // Mouse Position
  ImGui::Separator();
  ImGui::Text("Mouse Position: (%f, %f, %f)", mx, my, mz);
  ImApp::Pixel color = get_color(mcell);
  ImVec4 fcolor(static_cast<float>(color.r()) / 255.f,
                static_cast<float>(color.g()) / 255.f,
                static_cast<float>(color.b()) / 255.f,
                static_cast<float>(color.a()) / 255.f);

  if (!mcell) {
    ImGui::Text("Cell for given position is not defined.");
  } else {
    ImGui::Text("Cell ID: %i", mcell->id());
    ImGui::Text("Cell Instance: %i", mcell_inst);
    ImGui::Text("Cell Name: %s", mcell->name().data());

    if (colorby == ColorBy::Cell &&
        ImGui::ColorEdit3("", reinterpret_cast<float*>(&fcolor))) {
      color.r() = static_cast<uint8_t>(fcolor.x * 255.f);
      color.g() = static_cast<uint8_t>(fcolor.y * 255.f);
      color.b() = static_cast<uint8_t>(fcolor.z * 255.f);

      cell_id_to_color[mcell->id()] = color;

      must_rerender = true;
    }
  }
  if (!mmaterial) {
    ImGui::Text("Material for given position is not defined.");
  } else {
    ImGui::Text("Material ID: %i", mmaterial->id());
    ImGui::Text("Material Name: %s", mmaterial->name().data());

    if (colorby == ColorBy::Material &&
        ImGui::ColorEdit3("", reinterpret_cast<float*>(&fcolor))) {
      color.r() = static_cast<uint8_t>(fcolor.x * 255.f);
      color.g() = static_cast<uint8_t>(fcolor.y * 255.f);
      color.b() = static_cast<uint8_t>(fcolor.z * 255.f);

      material_id_to_color[mmaterial->id()] = color;

      must_rerender = true;
    }
  }

  // Save Image
  ImGui::Separator();
  if (ImGui::Button("Save Plot")) ImGui::OpenPopup("Save Plot");
  if (ImGui::BeginPopup("Save Plot")) {
    static char fn_str[500] = "";

    ImGui::InputTextWithHint("", "Enter File Name (i.e. reactor)", fn_str,
                             IM_ARRAYSIZE(fn_str));

    if (ImGui::Button("Save JPG")) {
      std::filesystem::path fname(fn_str);
      fname += ".jpg";
      image.save_jpg(fname);
      ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Save PNG")) {
      std::filesystem::path fname(fn_str);
      fname += ".png";
      image.save_png(fname);
      ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
  }

  ImGui::End();
}

void GuiPlotter::render_image() {
  {
    // Get the tracking direction
    const Direction u = this->get_tracking_direction();

// Go through each pixel
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (uint32_t i = 0; i < image.height(); i++) {
      Position strt = get_start_position(i);

      // Initalize the tracker
      Tracker trkr(strt, u);
      ::Cell* cell = trkr.cell();
      ImApp::Pixel pixel_color = this->get_color(cell);

      uint32_t j = 0;
      while (j < image.width()) {
        // Get the boundary
        auto bound = trkr.get_nearest_boundary();

        // Get the number of pixels till the boundary
        const double pixels_to_bound = bound.distance / dist_per_pixel;
        uint32_t npixels = static_cast<uint32_t>(std::round(pixels_to_bound));
        if (npixels > (image.width() - j) || bound.distance == INF) {
          npixels = image.width() - j;
        }
        const double npixels_dist =
            static_cast<double>(npixels) * dist_per_pixel;

        // Set all pixels
        for (uint32_t p = 0; p < npixels; p++) {
          if (j >= image.width()) break;
          if (outline_boundaries && npixels > 0 &&
              (p == npixels - 1 || j == image.width() - 1)) {
            image.at(i, j) = ImApp::Pixel(0, 0, 0);
          } else {
            image.at(i, j) = pixel_color;
          }
          j++;
        }
        if (j >= image.width()) break;

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
          if (j >= image.width()) break;
          image.at(i, j) = pixel_color;
          j++;
        }
      }  // while j < plot_width_
    }    // For i which is parallel
  }

  if (outline_boundaries) {
    // Get the tracking direction
    const Direction u = this->get_comp_tracking_direction();

    // Go through each pixel
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (uint32_t j = 0; j < image.width(); j++) {
      Position strt = get_comp_start_position(j);

      // Initalize the tracker
      Tracker trkr(strt, u);

      uint32_t i = 0;
      while (i < image.height()) {
        // Get the boundary
        auto bound = trkr.get_nearest_boundary();

        // Get the number of pixels till the boundary
        const double pixels_to_bound = bound.distance / dist_per_pixel;
        uint32_t npixels = static_cast<uint32_t>(std::round(pixels_to_bound));
        if (npixels > (image.height() - i) || bound.distance == INF) {
          npixels = image.height() - i;
        }
        const double npixels_dist =
            static_cast<double>(npixels) * dist_per_pixel;

        // Set all pixels
        for (uint32_t p = 0; p < npixels; p++) {
          if (i >= image.height()) break;
          if (npixels > 0 && (p == npixels - 1 || i == image.height() - 1)) {
            image.at(i, j) = ImApp::Pixel(0, 0, 0);
          }
          i++;
        }
        if (i >= image.height()) break;

        // Cross boundary, update cell, and get new pixel
        if (pixels_to_bound - static_cast<double>(npixels) < 0.5) {
          trkr.move(npixels_dist + dist_per_pixel);
        } else {
          trkr.move(npixels_dist);
        }
        trkr.restart_get_current();
        if (pixels_to_bound - static_cast<double>(npixels) < 0.5) {
          if (i >= image.height()) break;
          i++;
        }
      }  // while i < image.height
    }    // For j which is parallel
  }
}

ImApp::Pixel GuiPlotter::get_random_color() {
  uint8_t r = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
  uint8_t g = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));
  uint8_t b = static_cast<uint8_t>(255.0 * RNG::rand(settings::rng));

  return ImApp::Pixel(r, g, b);
}

ImApp::Pixel GuiPlotter::get_color(::Cell* cell) {
  ImApp::Pixel pixel_color = background;
  // If pointer isn't nullpntr, get pixel
  if (cell != nullptr) {
    if (colorby == ColorBy::Cell) {
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

Direction GuiPlotter::get_tracking_direction() const {
  switch (basis) {
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

Direction GuiPlotter::get_comp_tracking_direction() const {
  switch (basis) {
    case Basis::XY:
      return Direction(1., 0., 0.);
      break;
    case Basis::YZ:
      return Direction(0., 1., 0.);
      break;
    case Basis::XZ:
      return Direction(1., 0., 0.);
      break;
  }

  // NEVER GETS HERE
  return Direction(1., 0., 0.);
}

Position GuiPlotter::get_start_position(uint64_t i) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= image.height()) {
    fatal_error("Trying to deffine invalid pixel for plot.");
  }

  if (basis == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    const double dx = height / static_cast<double>(image.height());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double y_low = oy - (0.5 * width);

    // Get coordinate of pixel
    const double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    const double y = y_low;

    // Return coordinate
    return {x, y, oz};
  } else if (basis == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    const double dy = height / static_cast<double>(image.height());

    // Get upper corner off plot
    const double y_low = oy - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double y = (static_cast<double>(i) + 0.5) * dy + y_low;
    const double z = z_low;

    // Return coordinate
    return {ox, y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    const double dx = height / static_cast<double>(image.height());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    const double z = z_low;

    // Return coordinate
    return {x, oy, z};
  }
}

Position GuiPlotter::get_comp_start_position(uint64_t j) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (j >= image.width()) {
    fatal_error("Trying to deffine invalid pixel for plot.");
  }

  if (basis == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    const double dy = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double y_low = oy - (0.5 * width);

    // Get coordinate of pixel
    const double x = x_low;
    const double y = (static_cast<double>(j) + 0.5) * dy + y_low;

    // Return coordinate
    return {x, y, oz};
  } else if (basis == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    const double dz = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double y_low = oy - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double y = y_low;
    const double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {ox, y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    const double dz = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double x = x_low;
    const double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {x, oy, z};
  }
}

Position GuiPlotter::get_pixel_position(uint64_t i, uint64_t j) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= image.height() || j >= image.width()) {
    fatal_error("Trying to deffine invalid pixel for plot.");
  }

  if (basis == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    const double dx = height / static_cast<double>(image.height());
    const double dy = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double y_low = oy - (0.5 * width);

    // Get coordinate of pixel
    const double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    const double y = (static_cast<double>(j) + 0.5) * dy + y_low;

    // Return coordinate
    return {x, y, oz};
  } else if (basis == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    const double dy = height / static_cast<double>(image.height());
    const double dz = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double y_low = oy - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double y = (static_cast<double>(i) + 0.5) * dy + y_low;
    const double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {ox, y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    const double dx = height / static_cast<double>(image.height());
    const double dz = width / static_cast<double>(image.width());

    // Get upper corner off plot
    const double x_low = ox - (0.5 * height);
    const double z_low = oz - (0.5 * width);

    // Get coordinate of pixel
    const double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    const double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {x, oy, z};
  }
}

}  // namespace plotter

#endif  // ABEILLE_GUI_PLOT
