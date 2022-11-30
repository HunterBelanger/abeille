#ifndef GUI_PLOTTER_H
#define GUI_PLOTTER_H
#ifdef ABEILLE_GUI_PLOT

#include <ImApp/imapp.hpp>
#include <geometry/cell.hpp>
#include <map>
#include <mutex>
#include <utils/direction.hpp>
#include <utils/position.hpp>

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
