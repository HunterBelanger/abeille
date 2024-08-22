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
#ifndef CYLINDRICAL_OSCILLATION_NOISE_SOURCE_H
#define CYLINDRICAL_OSCILLATION_NOISE_SOURCE_H

#include <noise_source/oscillation_noise_source.hpp>
#include <utils/position.hpp>
#include <utils/tabulated_1d.hpp>

#include <yaml-cpp/yaml.h>

class CylindricalOscillationNoiseSource : public OscillationNoiseSource {
 public:
  CylindricalOscillationNoiseSource(Position origin, double len, double rad,
                                    char ax, const pndl::Tabulated1D& eps_tot,
                                    const pndl::Tabulated1D& eps_fis,
                                    const pndl::Tabulated1D& eps_sct,
                                    double angular_frequency, double phase);

  bool is_inside(const Position& r) const override final;
  std::complex<double> dEt(const Position& r, double E,
                           double w) const override final;

  std::complex<double> dEt_Et(const Position& r, double E,
                              double w) const override final;
  std::complex<double> dEf_Ef(const Position& r, double E,
                              double w) const override final;
  std::complex<double> dEelastic_Eelastic(const Position& r, double E,
                                          double w) const override final;
  std::complex<double> dEmt_Emt(uint32_t mt, const Position& r, double E,
                                double w) const override final;

 private:
  Position origin_;
  double length_, radius_;
  char axis_;
  double w0_, phase_;
  pndl::Tabulated1D eps_t_;
  pndl::Tabulated1D eps_f_;
  pndl::Tabulated1D eps_s_;
};

std::shared_ptr<OscillationNoiseSource>
make_cylindrical_oscillation_noise_source(const YAML::Node& snode);

#endif
