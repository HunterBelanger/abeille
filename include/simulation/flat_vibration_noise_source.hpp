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
#ifndef FLAT_VIBRATION_NOISE_SOURCE_H
#define FLAT_VIBRATION_NOISE_SOURCE_H

#include <materials/material.hpp>
#include <simulation/vibration_noise_source.hpp>

#include <complex>
#include <memory>
#include <unordered_map>

class FlatVibrationNoiseSource : public VibrationNoiseSource {
 public:
  enum class Basis { X, Y, Z };

 public:
  FlatVibrationNoiseSource(Position low, Position hi, Basis basis,
                           std::shared_ptr<Material> mat_pos,
                           std::shared_ptr<Material> mat_neg,
                           double angular_frequency);

  bool is_inside(const Position& r) const override final;

  std::complex<double> dEt(const Position& r, double E,
                           double w) const override final;

  std::complex<double> dEt_Et(const Position& r, double E,
                              double w) const override final;

  std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                          double w) const override final;

 private:
  Position low_, hi_;
  Basis basis_;
  std::shared_ptr<Material> material_pos_;  // Material on the positive side
  std::shared_ptr<Material> material_neg_;  // Material on the negative side
  double x0_;                               // Midpoint of interface
  double w0_;                               // Angular frequency of vibration
  double eps_;                              // Magnitude of oscillation
  std::unordered_map<uint32_t, double> Delta_N;  // Contians N_neg - N_pos

  // Function to get xs
  double Et(double x, double E) const;
  // Function to get Ei_neg(E) - Ei_pos(E)
  double Delta_Et(double E) const;

  std::complex<double> C_R(uint32_t n, double x) const;
  std::complex<double> C_L(uint32_t n, double x) const;

  double get_x(const Position& r) const;

  bool negative_material(double x) const;

  double get_nuclide_concentration(const Material& mat,
                                   uint32_t nuclide_id) const;

  static constexpr std::complex<double> i{0., 1.};
};

std::shared_ptr<VibrationNoiseSource> make_flat_vibration_noise_source(
    const YAML::Node& snode);

#endif
