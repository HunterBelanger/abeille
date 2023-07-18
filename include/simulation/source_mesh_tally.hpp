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
#ifndef SOURCE_MESH_TALLY_H
#define SOURCE_MESH_TALLY_H

#include <simulation/mesh_tally.hpp>
#include <simulation/particle.hpp>

#include <yaml-cpp/yaml.h>

class SourceMeshTally : public MeshTally {
 public:
  enum class Quantity { Source, RealSource, ImagSource };

  SourceMeshTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                  uint64_t nz, const std::vector<double>& ebounds, Quantity q,
                  std::string fname)
      : MeshTally(low, hi, nx, ny, nz, ebounds, fname), quantity(q) {}

  void score_source(const BankedParticle& p);

  bool noise_like_score() const {
    if (quantity == Quantity::RealSource || quantity == Quantity::ImagSource)
      return true;
    return false;
  }

  std::string estimator_str() const override final { return "source"; }

  std::string quantity_str() const override final;

  std::uint32_t mt() const override final { return 0; }

 private:
  Quantity quantity;
};

std::shared_ptr<SourceMeshTally> make_source_mesh_tally(const YAML::Node& node);

#endif
