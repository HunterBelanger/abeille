#ifdef ABEILLE_GUI_PLOT

#include <plotting/gui_plotter.hpp>

namespace plotter {

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

  this->render_controls();
  this->render_viewport();
}

void GuiPlotter::render_viewport() {
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Viewport");
  ImGui::End();
}

void GuiPlotter::render_controls() {
  ImGui::SetNextWindowSize({500, 500}, ImGuiCond_Once);
  ImGui::Begin("Controls");
  ImGui::End();
}
  
}

#endif // ABEILLE_GUI_PLOT
