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
#ifndef NOISE_SOURCE_H
#define NOISE_SOURCE_H

#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

#include <complex>

// Pure virtual interface to allow for the sampling of neutron noise particles.
class NoiseSource {
 public:
  NoiseSource() = default;
  virtual ~NoiseSource() = default;

  virtual bool is_inside(const Position& r) const = 0;

  virtual std::complex<double> dEt(const Position& r, double E,
                                   double w) const = 0;

  virtual std::complex<double> dEt_Et(const Position& r, double E,
                                      double w) const = 0;
};

#endif
