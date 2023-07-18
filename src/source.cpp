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
#include <simulation/source.hpp>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

Source::Source(std::shared_ptr<SpatialDistribution> spatial,
               std::shared_ptr<DirectionDistribution> direction,
               std::shared_ptr<EnergyDistribution> energy, bool fissile_only,
               double weight)
    : spatial_(spatial),
      direction_(direction),
      energy_(energy),
      fissile_only_(fissile_only),
      weight_(weight) {
  if (weight_ <= 0.) {
    fatal_error("Source weight must be greater than zero.");
  }
}

Particle Source::generate_particle(pcg32& rng) const {
  // Sample Direction
  Direction u = direction_->sample(rng);

  // Sample energy
  double E = energy_->sample(rng);
  int E_count = 0;
  do {
    if (E_count > 200) {
      fatal_error("Exceded 200 samplings of energy.");
    }

    E = energy_->sample(rng);
    E_count++;
  } while (E <= settings::min_energy || settings::max_energy <= E);

  // Sample position
  Position r = spatial_->sample(rng);
  Tracker trkr(r, u);

  // Keep sampling new position untill a position inside the geometry
  // is found.
  while (trkr.is_lost()) {
    r = spatial_->sample(rng);
    trkr.set_r(r);
    trkr.restart_get_current();
  }

  // If we only want positions inside of fissile regions, then we need
  // to check this.
  if (fissile_only_) {
    int count = 0;
    while (trkr.material()->fissile() == false || trkr.is_lost()) {
      if (count == 201) {
        fatal_error(
            "Exceded 200 samplings of position in fissile-only source.");
      }

      r = spatial_->sample(rng);
      trkr.set_r(r);
      trkr.restart_get_current();
      count++;
    }
  }

  return Particle(r, u, E, 1.0);
}

std::shared_ptr<Source> make_source(const YAML::Node& source_node) {
  // Make sure it is a map
  if (!source_node.IsMap()) {
    fatal_error("Source entry must be a map.");
  }

  // Get weight
  if (!source_node["weight"] || !source_node["weight"].IsScalar()) {
    fatal_error("No weight given to source.");
  }
  double wgt = source_node["weight"].as<double>();

  // Get fissile only
  bool fissile_only = false;
  if (source_node["fissile-only"] && source_node["fissile-only"].IsScalar()) {
    fissile_only = source_node["fissile-only"].as<bool>();
  } else if (source_node["fissile-only"]) {
    fatal_error("Invalid fissile-only entry in source description.");
  }

  // Get spatial distribution
  if (!source_node["spatial"] || !source_node["spatial"].IsMap()) {
    fatal_error(
        "No valid spatial distribution entry provided for source "
        "distribution.");
  }
  std::shared_ptr<SpatialDistribution> spatial =
      make_spatial_distribution(source_node["spatial"]);

  // Get direction distribution
  if (!source_node["direction"] || !source_node["direction"].IsMap()) {
    fatal_error(
        "No valid direction distribution entry provided for source "
        "distribution.");
  }
  std::shared_ptr<DirectionDistribution> direction =
      make_direction_distribution(source_node["direction"]);

  // Get energy distribution
  if (!source_node["energy"] || !source_node["energy"].IsMap()) {
    fatal_error(
        "No valid energy distribution entry provided for source distribution.");
  }
  std::shared_ptr<EnergyDistribution> energy =
      make_energy_distribution(source_node["energy"]);

  return std::make_shared<Source>(spatial, direction, energy, fissile_only,
                                  wgt);
}
