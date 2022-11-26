#ifndef GUI_PLOTTER_H
#define GUI_PLOTTER_H
#ifdef ABEILLE_GUI_PLOT

#include <ImApp/imapp.hpp>

namespace plotter {
  
  class GuiPlotter : public ImApp::Layer {
    public:
      GuiPlotter() = default;

      void render() override final;

    private:
      void render_viewport();
      void render_controls();
  };

}

#endif // ABEILLE_GUI_PLOT
#endif // GUI_PLOTTER_H
