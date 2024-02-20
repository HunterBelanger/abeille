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
#ifndef SOURCE_H
#define SOURCE_H

#include <simulation/particle.hpp>
#include <source/direction_distribution.hpp>
#include <source/energy_distribution.hpp>
#include <source/spatial_distribution.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class Source {
 public:
  Source(std::shared_ptr<SpatialDistribution> spatial,
         std::shared_ptr<DirectionDistribution> direction,
         std::shared_ptr<EnergyDistribution> energy, bool fissile_only,
         double weight);

  Particle generate_particle(RNG& rng) const;
  double wgt() const { return weight_; }
  bool fissile_only() const { return fissile_only_; }

 private:
  std::shared_ptr<SpatialDistribution> spatial_;
  std::shared_ptr<DirectionDistribution> direction_;
  std::shared_ptr<EnergyDistribution> energy_;
  bool fissile_only_;
  double weight_;
};

// Helper function to build a source entry
std::shared_ptr<Source> make_source(const YAML::Node& source_node);

#endif
