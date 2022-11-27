#ifdef ABEILLE_GUI_PLOT

#include <plotting/gui_plotter.hpp>
#include <plotting/pixel.hpp>
#include <utils/error.hpp>
#include <simulation/tracker.hpp>

namespace plotter {

// Maps for colors from plotter.hpp
extern std::map<uint32_t, Pixel> cell_id_to_color;
extern std::map<uint32_t, Pixel> material_id_to_color;

GuiPlotter::GuiPlotter():
 cell_id_to_color(),
 material_id_to_color(),
 image(500, 500),
 create_color_mutex(),
 height(10.),
 width(10.),
 ox(0.), oy(0.), oz(0.),
 background(),
 rng(),
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
  if (ImGui::BeginMainMenuBar()) {
    if (ImGui::BeginMenu("Options")) {
      if (ImGui::BeginMenu("Save")) {
        ImGui::EndMenu();
      }
      ImGui::EndMenu();
    }
    ImGui::EndMainMenuBar();
  }

  // We put a Dockspace over the entire viewport.
  ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

  this->render_viewport();
  this->render_controls();
}

void GuiPlotter::render_viewport() {
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Viewport");
  
  // First, get the window size. If the size doesn't match the current
  // image size, we must set must_rerender to true. 
  std::uint32_t wwidth = static_cast<std::uint32_t>(ImGui::GetWindowWidth());
  std::uint32_t wheight = static_cast<std::uint32_t>(ImGui::GetWindowHeight());
  if ((wwidth != image.width()) ||
      (wheight != image.height())) {
    image.resize(wheight, wwidth);

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

  // Now we need to add the image to the window
  void* texture_id = reinterpret_cast<void*>(
      static_cast<intptr_t>(image.ogl_texture_id().value()));
  ImGui::Image(texture_id, ImVec2(static_cast<float>(image.width()),
        static_cast<float>(image.height())));

  ImGui::End();
}

void GuiPlotter::render_controls() {
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Controls");

  // Origin
  ImGui::Separator();
  ImGui::Text("Plot Origin");
  if(ImGui::InputDouble("X [cm] :", &ox, 0., 0., "%E")) must_rerender = true;
  if(ImGui::InputDouble("Y [cm] :", &oy, 0., 0., "%E")) must_rerender = true;
  if(ImGui::InputDouble("Z [cm] :", &oz, 0., 0., "%E")) must_rerender = true;

  // Physical Dimensions of plot
  ImGui::Separator();
  ImGui::Text("Width/Height");
  ImGui::RadioButton("Width", &adjust_w_or_h, 0); ImGui::SameLine();
  ImGui::RadioButton("Height", &adjust_w_or_h, 1);
  if (adjust_w_or_h == 0) {
    if (ImGui::InputDouble("Width [cm] :", &width, 1.E-6, 1.E32, "%E")) {
      must_rerender = true;

      // Must recalculate height
      dist_per_pixel = width / static_cast<double>(image.width());
      height = static_cast<double>(image.height()) * dist_per_pixel;
    }
  } else if(adjust_w_or_h == 1) {
    if(ImGui::InputDouble("Height [cm] :", &height, 1.E-6, 1.E32, "%E")) {
      must_rerender = true; 
      
      // Must recalculate width
      dist_per_pixel = height / static_cast<double>(image.height());
      width = static_cast<double>(image.width()) * dist_per_pixel;
    }
  }
  ImGui::Text("Width [cm]: %E, Height [cm]: %E", width, height);

  // Slice Basis
  ImGui::Separator();
  ImGui::Text("Slice Basis");
  if(ImGui::RadioButton("XY", reinterpret_cast<int*>(&basis), Basis::XY)) must_rerender = true;
  ImGui::SameLine();
  if(ImGui::RadioButton("XZ", reinterpret_cast<int*>(&basis), Basis::XZ)) must_rerender = true;
  ImGui::SameLine();
  if(ImGui::RadioButton("YZ", reinterpret_cast<int*>(&basis), Basis::YZ)) must_rerender = true;

  // Color Method
  ImGui::Separator();
  ImGui::Text("Coloring");
  if(ImGui::RadioButton("Cell", reinterpret_cast<int*>(&colorby), ColorBy::Cell)) must_rerender = true;
  ImGui::SameLine();
  if(ImGui::RadioButton("Material", reinterpret_cast<int*>(&colorby), ColorBy::Material)) must_rerender = true;

  if (ImGui::Checkbox("Mark Boundaries", &outline_boundaries)) must_rerender = true;

  ImGui::End();
}

void GuiPlotter::render_image() {
  // Get the tracking direction
  const Direction u = this->get_tracking_direction();

// Go through each pixel
#ifdef _OPENMP
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
      const double npixels_dist = static_cast<double>(npixels) * dist_per_pixel;

      // Set all pixels
      for (uint32_t p = 0; p < npixels; p++) {
        if (j >= image.width()) break;
        if (outline_boundaries && npixels > 0 &&
            (p == npixels-1 || j == image.width()-1)) {
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

ImApp::Pixel GuiPlotter::get_random_color() {
  uint8_t r = static_cast<uint8_t>(255.0 * RNG::rand(rng));
  uint8_t g = static_cast<uint8_t>(255.0 * RNG::rand(rng));
  uint8_t b = static_cast<uint8_t>(255.0 * RNG::rand(rng));

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

Position GuiPlotter::get_start_position(uint64_t i) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= image.height()) {
    std::string mssg = "Trying to deffine invalid pixel for plot.";
    fatal_error(mssg, __FILE__, __LINE__);
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

Position GuiPlotter::get_pixel_position(uint64_t i, uint64_t j) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= image.height() || j >= image.width()) {
    std::string mssg = "Trying to deffine invalid pixel for plot.";
    fatal_error(mssg, __FILE__, __LINE__);
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
  
}

#endif // ABEILLE_GUI_PLOT
