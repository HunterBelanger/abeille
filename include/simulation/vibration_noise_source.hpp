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
#ifndef VIBRATION_NOISE_SOURCE_H
#define VIBRATION_NOISE_SOURCE_H

#include <simulation/noise_source.hpp>

#include <complex>
#include <memory>
#include <unordered_map>
#include <vector>

// Pure virtual interface to allow for the sampling of neutron noise particles.
class VibrationNoiseSource : public NoiseSource {
 public:
  struct NuclideInfo {
    uint32_t id;
    double concentration;
  };

  VibrationNoiseSource() : nuclides_(), nuclide_info_() {}
  virtual ~VibrationNoiseSource() = default;

  virtual std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                                  double w) const = 0;

  const std::vector<uint32_t>& nuclides() const { return nuclides_; }

  const std::unordered_map<uint32_t, NuclideInfo>& nuclide_info() const {
    return nuclide_info_;
  }

 protected:
  std::vector<uint32_t> nuclides_;
  std::unordered_map<uint32_t, NuclideInfo> nuclide_info_;
};

#endif
