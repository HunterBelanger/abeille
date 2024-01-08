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
#include <cancelator/approximate_mesh_cancelator.hpp>
#include <cancelator/cancelator.hpp>
#include <cancelator/exact_mg_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

std::shared_ptr<Cancelator> make_cancelator(const YAML::Node& node) {
  if (!node["type"] || !node["type"].IsScalar()) {
    fatal_error("Invalid type entry for cancelator.");
  }
  std::string type = node["type"].as<std::string>();

  std::shared_ptr<Cancelator> cancelator = nullptr;

  if (type == "approximate") {
    cancelator = make_approximate_mesh_cancelator(node);
  } else if (type == "exact") {
    // Check that we are not in CE mode !
    if (settings::energy_mode == settings::EnergyMode::CE) {
      fatal_error(
          "exact cancelators are not currently supported for continuous "
          "energy.");
    }

    // looks like we meet all requirments, make the cancelator
    cancelator = make_exact_mg_cancelator(node);
  } else {
    fatal_error("Unkown cancelator type \"" + type + "\".");
  }

  return cancelator;
}
