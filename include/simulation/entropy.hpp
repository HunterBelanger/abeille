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
#ifndef ENTROPY_H
#define ENTROPY_H

#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

class Entropy {
 public:
  enum Sign { Positive, Negative, Total };

 public:
  Entropy(Position low, Position up, std::array<uint32_t, 3> shp,
          Sign sgn = Sign::Positive)
      : lower_corner{low},
        upper_corner{up},
        shape{shp},
        nbins{0},
        dx{},
        dy{},
        dz{},
        bins{},
        total_weight{},
        sign_{sgn} {
    nbins = shape[0] * shape[1] * shape[2];
    bins.resize(nbins, 0.);
    total_weight = 0.;
    dx = (upper_corner.x() - lower_corner.x()) / static_cast<double>(shape[0]);
    dy = (upper_corner.y() - lower_corner.y()) / static_cast<double>(shape[1]);
    dz = (upper_corner.z() - lower_corner.z()) / static_cast<double>(shape[2]);
  }
  ~Entropy() = default;

  Sign sign() const { return sign_; }

  void set_sign(Sign sgn) { sign_ = sgn; }

  void add_point(const Position& r, const double& w);

  void synchronize_entropy_across_nodes();

  double calculate_entropy() const;

  double calculate_empty_fraction() const;

  void zero();

  double total_wgt() const { return total_weight; }

 private:
  Position lower_corner;
  Position upper_corner;
  std::array<std::uint32_t, 3> shape;
  std::uint32_t nbins;
  double dx, dy, dz;
  std::vector<double> bins;
  double total_weight;
  Sign sign_;

};  // Entropy

Entropy make_entropy(const YAML::Node& entropy);

#endif  // MG_ENTROPY_H
