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
#include <source/mono_energetic.hpp>
#include <utils/error.hpp>

#include <memory>

MonoEnergetic::MonoEnergetic(double energy) : energy_(energy) {
  if (energy_ <= 0.) {
    fatal_error("Cannot have a negative or zero energy.");
  }
}

double MonoEnergetic::sample(pcg32& /*rng*/) const { return energy_; }

std::shared_ptr<MonoEnergetic> make_mono_energetic_distribution(
    const YAML::Node& node) {
  if (!node["energy"] || !node["energy"].IsScalar()) {
    fatal_error("No valid energy entry to mono-energetic distribution.");
  }

  double E = node["energy"].as<double>();

  return std::make_shared<MonoEnergetic>(E);
}
