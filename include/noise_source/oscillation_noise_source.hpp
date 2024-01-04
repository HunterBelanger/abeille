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
#ifndef OSCILLATION_NOISE_SOURCE_H
#define OSCILLATION_NOISE_SOURCE_H

#include <noise_source/noise_source.hpp>
#include <utils/position.hpp>

// Pure virtual interface to allow for the sampling of neutron noise particles.
class OscillationNoiseSource : public NoiseSource {
 public:
  OscillationNoiseSource() = default;
  virtual ~OscillationNoiseSource() = default;

  virtual std::complex<double> dEt_Et(const Position& r, double E,
                                      double w) const = 0;

  virtual std::complex<double> dEf_Ef(const Position& r, double E,
                                      double w) const = 0;

  virtual std::complex<double> dEelastic_Eelastic(const Position& r, double E,
                                                  double w) const = 0;

  virtual std::complex<double> dEmt_Emt(uint32_t mt, const Position& r,
                                        double E, double w) const = 0;
};

#endif
