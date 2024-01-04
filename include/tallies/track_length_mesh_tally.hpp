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
#ifndef TRACK_LENGTH_MESH_TALLY_H
#define TRACK_LENGTH_MESH_TALLY_H

#include <tallies/mesh_tally.hpp>

#include <yaml-cpp/yaml.h>

class TrackLengthMeshTally : public MeshTally {
 public:
  TrackLengthMeshTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                       uint64_t nz, const std::vector<double>& ebounds,
                       Quantity q, std::string fname, uint32_t mt = 0)
      : MeshTally(low, hi, nx, ny, nz, ebounds, fname), quantity(q), mt_(mt) {}

  void score_flight(const Particle& p, double d, MaterialHelper& mat);

  std::string estimator_str() const override final { return "track-length"; }

  std::string quantity_str() const override final;

  std::uint32_t mt() const override final { return mt_; }

 private:
  Quantity quantity;
  uint32_t mt_;

  void initialize_indices(const Position& r, const Direction& u, int& i, int& j,
                          int& k, std::array<int, 3>& on);
  std::pair<double, int> distance_to_next_index(const Position& r,
                                                const Direction& u, int i,
                                                int j, int k,
                                                const std::array<int, 3>& on);
  bool find_entry_point(Position& r, const Direction& u,
                        double& d_flight) const;
  void update_indices(int key, int& i, int& j, int& k, std::array<int, 3>& on);

  double get_base_score(const Particle& p, MaterialHelper& mat) const;
};

std::shared_ptr<TrackLengthMeshTally> make_track_length_mesh_tally(
    const YAML::Node& node);

#endif
