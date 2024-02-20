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
#ifndef CONE_H
#define CONE_H

#include <source/direction_distribution.hpp>

#include <yaml-cpp/yaml.h>

//         u
//         ^
//  \      |  ap  /
//   \     |-----/
//    \    |    /
//     \   |   /
//      \  |  /
//       \ | /
//        \|/
//         *

class Cone : public DirectionDistribution {
 public:
  // aperture given in radians (0 is a MonoDirectional beam). It is the aperture
  // of HALF the beam.
  Cone(Direction u, double aperture);

  Direction sample(RNG& rng) const override final;

  Direction direction() const { return direction_; }

  // Returns aperture in radians.
  double aperture() const;

 private:
  Direction direction_;
  double aperture_;
};

std::shared_ptr<Cone> make_cone_distribution(const YAML::Node& node);

#endif
