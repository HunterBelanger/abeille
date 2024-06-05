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
#include <materials/material_helper.hpp>
#include <noise_source/square_oscillation_noise_source.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/direction.hpp>
#include <utils/error.hpp>

#include <sstream>

SquareOscillationNoiseSource::SquareOscillationNoiseSource(
    Position low, Position hi, double eps_tot, double eps_fis, double eps_sct,
    double angular_frequency, double phase)
    : low_(low),
      hi_(hi),
      w0_(angular_frequency),
      phase_(phase),
      eps_t_(eps_tot),
      eps_f_(eps_fis),
      eps_s_(eps_sct) {
  // Check low and high
  if (low_.x() >= hi_.x() || low_.y() >= hi_.y() || low_.z() >= hi_.z()) {
    fatal_error(
        "Low is greater than or equal to hi in SquareOscillationNoiseSource.");
  }

  // Make sure frequency is positive
  if (w0_ <= 0.) {
    fatal_error("Negative or zero noise frequency.");
  }

  // Make sure epsilons are positive
  if (eps_t_ <= 0.) {
    fatal_error("Negative or zero epsilon total.");
  }

  if (eps_f_ <= 0.) {
    fatal_error("Negative or zero epsilon fission.");
  }

  if (eps_s_ <= 0.) {
    fatal_error("Negative or zero epsilon scatter.");
  }
}

bool SquareOscillationNoiseSource::is_inside(const Position& r) const {
  if (r.x() > low_.x() && r.y() > low_.y() && r.z() > low_.z() &&
      r.x() < hi_.x() && r.y() < hi_.y() && r.z() < hi_.z())
    return true;

  return false;
}

std::complex<double> SquareOscillationNoiseSource::dEt(const Position& r,
                                                       double E,
                                                       double w) const {
  // Get the material
  Tracker trkr(r, {1., 0., 0.});
  auto material = trkr.material();

  // Make sure material is valid
  if (!material) {
    std::stringstream mssg;
    mssg << "No material found at position " << r << ".";
    fatal_error(mssg.str());
  }

  // Make material helper and evaluate XS
  MaterialHelper mat(material, E);
  double xs = mat.Et(E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * xs * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * xs * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> SquareOscillationNoiseSource::dEt_Et(const Position& /*r*/,
                                                          double /*E*/,
                                                          double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> SquareOscillationNoiseSource::dEf_Ef(const Position& /*r*/,
                                                          double /*E*/,
                                                          double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_ * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_ * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> SquareOscillationNoiseSource::dEelastic_Eelastic(
    const Position& /*r*/, double /*E*/, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> SquareOscillationNoiseSource::dEmt_Emt(
    uint32_t /*mt*/, const Position& /*r*/, double /*E*/, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::shared_ptr<OscillationNoiseSource> make_square_oscillation_noise_source(
    const YAML::Node& snode) {
  // Get low
  if (!snode["low"] || !snode["low"].IsSequence() ||
      !(snode["low"].size() == 3)) {
    fatal_error("No valid low entry for square oscillation noise source.");
  }

  double xl = snode["low"][0].as<double>();
  double yl = snode["low"][1].as<double>();
  double zl = snode["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!snode["hi"] || !snode["hi"].IsSequence() || !(snode["hi"].size() == 3)) {
    fatal_error("No valid hi entry for square oscillation noise source.");
  }

  double xh = snode["hi"][0].as<double>();
  double yh = snode["hi"][1].as<double>();
  double zh = snode["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get frequency
  if (!snode["angular-frequency"] || !snode["angular-frequency"].IsScalar()) {
    fatal_error(
        "No valid angular-frequency entry for square oscillation noise "
        "source.");
  }

  double w0 = snode["angular-frequency"].as<double>();

  if (w0 <= 0.) {
    fatal_error(
        "Angular frequency can not be negative or zero for square oscillation "
        "noise source.");
  }

  // Get phase (is present)
  double phase = 0.;
  if (snode["phase"] && snode["phase"].IsScalar()) {
    phase = snode["phase"].as<double>();
  } else if (snode["phase"]) {
    fatal_error("Invalid phase entry for square oscillation noise source.");
  }

  // Get epsilons
  if (!snode["epsilon-total"] || !snode["epsilon-total"].IsScalar()) {
    fatal_error(
        "No valid epsilon-total entry for square oscillation noise source.");
  }

  double eps_t = snode["epsilon-total"].as<double>();

  if (eps_t <= 0.) {
    fatal_error(
        "Epsilon total can not be negative or zero for square oscillation "
        "noise source.");
  }

  if (!snode["epsilon-fission"] || !snode["epsilon-fission"].IsScalar()) {
    fatal_error(
        "No valid epsilon-fission entry for square oscillation noise source.");
  }

  double eps_f = snode["epsilon-fission"].as<double>();

  if (eps_f <= 0.) {
    fatal_error(
        "Epsilon fission can not be negative or zero for square oscillation "
        "noise source.");
  }

  if (!snode["epsilon-scatter"] || !snode["epsilon-scatter"].IsScalar()) {
    fatal_error(
        "No valid epsilon-scatter entry for square oscillation noise source.");
  }

  double eps_s = snode["epsilon-scatter"].as<double>();

  if (eps_s <= 0.) {
    fatal_error(
        "Epsilon scatter can not be negative or zero for square oscillation "
        "noise source.");
  }

  return std::make_shared<SquareOscillationNoiseSource>(
      r_low, r_hi, eps_t, eps_f, eps_s, w0, phase);
}
