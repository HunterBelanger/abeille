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
#ifndef SPATIAL_DISTRIBUTION_H
#define SPATIAL_DISTRIBUTION_H

#include <utils/position.hpp>
#include <utils/rng.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class SpatialDistribution {
 public:
  SpatialDistribution() {}
  virtual ~SpatialDistribution() = default;

  virtual Position sample(RNG& rng) const = 0;
};

// Helper function to build spatial distribution
std::shared_ptr<SpatialDistribution> make_spatial_distribution(
    const YAML::Node& node);

#endif
