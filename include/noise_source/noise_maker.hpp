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
#ifndef NOISE_MAKER_H
#define NOISE_MAKER_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <noise_source/oscillation_noise_source.hpp>
#include <noise_source/vibration_noise_source.hpp>
#include <utils/noise_parameters.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <vector>

class NoiseMaker {
 public:
  NoiseMaker() {}

  void add_noise_source(const YAML::Node& snode);

  void add_noise_source(std::shared_ptr<VibrationNoiseSource> ns) {
    vibration_noise_sources_.push_back(ns);
  }

  void add_noise_source(std::shared_ptr<OscillationNoiseSource> ns) {
    oscillation_noise_sources_.push_back(ns);
  }

  std::size_t num_noise_sources() const {
    return vibration_noise_sources_.size() + oscillation_noise_sources_.size();
  }

  NoiseParameters& noise_parameters() { return noise_params_; }

  const NoiseParameters& noise_parameters() const { return noise_params_; }

  void sample_noise_source(Particle& p, MaterialHelper& mat) const;

 private:
  std::vector<std::shared_ptr<VibrationNoiseSource>> vibration_noise_sources_;
  std::vector<std::shared_ptr<OscillationNoiseSource>>
      oscillation_noise_sources_;
  NoiseParameters noise_params_;

  bool is_inside(const Particle& p) const;
  std::complex<double> dEt(const Particle& p, double w) const;
  std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                          double w) const;
  std::unique_ptr<Material> make_fake_material(const Particle& p) const;

  void sample_noise_copy(Particle& p, MaterialHelper& mat,
                         const double w) const;

  void sample_vibration_noise_source(Particle& p, MaterialHelper& mat,
                                     const double keff, const double w) const;
  void sample_vibration_noise_fission(Particle& p, const Nuclide& nuclide,
                                      const MicroXSs& microxs,
                                      const std::complex<double>& dN_N,
                                      const double Etfake_Et, const double keff,
                                      const double w) const;
  void sample_vibration_noise_scatter(Particle& p, const Nuclide& nuclide,
                                      const MicroXSs& microxs,
                                      const std::complex<double>& dN_N,
                                      const double Etfake_Et,
                                      const double P_scatter) const;

  void sample_oscillation_noise_source(Particle& p, MaterialHelper& mat,
                                       const double keff, const double w) const;
  void sample_oscillation_noise_fission(Particle& p, const Nuclide& nuclide,
                                        const MicroXSs& microxs,
                                        const double keff,
                                        const double w) const;
  void sample_oscillation_noise_scatter(Particle& p, const Nuclide& nuclide,
                                        const MicroXSs& microxs,
                                        const double P_scatter,
                                        const double w) const;
};

#endif
