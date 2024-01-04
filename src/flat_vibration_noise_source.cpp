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
#include <noise_source/flat_vibration_noise_source.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>

FlatVibrationNoiseSource::FlatVibrationNoiseSource(
    Position low, Position hi, Basis basis, std::shared_ptr<Material> mat_pos,
    std::shared_ptr<Material> mat_neg, double angular_frequency)
    : low_(low),
      hi_(hi),
      basis_(basis),
      material_pos_(mat_pos),
      material_neg_(mat_neg),
      x0_(),
      w0_(angular_frequency),
      eps_(),
      Delta_N() {
  // Check low and high
  if (low_.x() >= hi_.x() || low_.y() >= hi_.y() || low_.z() >= hi_.z()) {
    fatal_error(
        "Low is greater than or equal to hi in FlatVibrationNoiseSource.");
  }

  // Make sure frequency is positive
  if (w0_ <= 0.) {
    fatal_error(
        "Negative or zero frequency provided to FlatVibrationNoiseSource.");
  }

  // Make sure materials aren't nullptr
  if (material_pos_ == nullptr) {
    fatal_error("Nullptr positive material given to FlatVibrationNoiseSource.");
  }

  if (material_neg_ == nullptr) {
    fatal_error("Nullptr negative material given to FlatVibrationNoiseSource.");
  }

  // Get x0_ and eps_, based on the basis
  switch (basis_) {
    case Basis::X:
      x0_ = 0.5 * (low_.x() + hi_.x());
      eps_ = (hi_.x() - low_.x()) / 2.;
      break;

    case Basis::Y:
      x0_ = 0.5 * (low_.y() + hi_.y());
      eps_ = (hi_.y() - low_.y()) / 2.;
      break;

    case Basis::Z:
      x0_ = 0.5 * (low_.z() + hi_.z());
      eps_ = (hi_.z() - low_.z()) / 2.;
      break;
  }

  // No we need to fill all of the nuclide info for NoiseSource
  // and for FlatVibrationNoiseSource. We start by filling out
  // nuclide_info_, and nuclides_.
  // Negative Material
  for (const auto& comp : material_neg_->components()) {
    auto nuc_id = comp.nuclide->id();
    nuclide_info_[nuc_id] = {nuc_id, comp.atoms_bcm};
    nuclides_.push_back(nuc_id);
  }
  // Positive Material
  for (const auto& comp : material_pos_->components()) {
    auto nuc_id = comp.nuclide->id();
    // Check if nuclide was already added
    if (nuclide_info_.find(nuc_id) != nuclide_info_.end()) {
      // Nuclide already added once. Update it's concentration
      // to contain the average concentration of the two.
      nuclide_info_[nuc_id].concentration += comp.atoms_bcm;
      nuclide_info_[nuc_id].concentration /= 2.;
    } else {
      // New nuclide
      nuclide_info_[nuc_id] = {nuc_id, comp.atoms_bcm};
      nuclides_.push_back(nuc_id);
    }
  }

  // Now we need to get Delta_N for each nuclide in the two materials.
  for (const auto& nuc : nuclides_) {
    const double N_neg = get_nuclide_concentration(*material_neg_, nuc);
    const double N_pos = get_nuclide_concentration(*material_pos_, nuc);
    Delta_N[nuc] = N_neg - N_pos;
  }

  // We now need to sort the nuclides list ! This is important for
  // efficiently getting the union of all nuclides in the NoiseMaker.
  std::sort(nuclides_.begin(), nuclides_.end());
}

double FlatVibrationNoiseSource::get_nuclide_concentration(
    const Material& mat, uint32_t nuclide_id) const {
  for (const auto& comp : mat.components()) {
    if (comp.nuclide->id() == nuclide_id) return comp.atoms_bcm;
  }

  return 0.;
}

double FlatVibrationNoiseSource::Et(double x, double E) const {
  if (x < x0_) {
    MaterialHelper mat_neg(material_neg_.get(), E);
    return mat_neg.Et(E);
  } else {
    MaterialHelper mat_pos(material_pos_.get(), E);
    return mat_pos.Et(E);
  }
}

double FlatVibrationNoiseSource::Delta_Et(double E) const {
  // Make material helpers and evaluate XS
  MaterialHelper mat_neg(material_neg_.get(), E);
  double xs_neg = mat_neg.Et(E);

  MaterialHelper mat_pos(material_pos_.get(), E);
  double xs_pos = mat_pos.Et(E);

  return xs_neg - xs_pos;
}

std::complex<double> FlatVibrationNoiseSource::C_R(uint32_t n, double x) const {
  const double dn = static_cast<double>(n);

  double rel_diff = (x - x0_) / eps_;
  if (rel_diff > 1.)
    rel_diff = 1.;
  else if (rel_diff < -1.)
    rel_diff = -1.;

  switch (n) {
    case 0:
      return {PI - 2. * std::asin(rel_diff), 0.};
      break;

    case 1:
      return {0., -2. * std::sqrt(1. - (rel_diff * rel_diff))};
      break;

    case 2:
      return {-2. * rel_diff * std::sqrt(1. - (rel_diff * rel_diff)), 0.};
      break;

    default:
      return (2. / dn) * std::sin(dn * std::acos(rel_diff)) *
             std::exp(-i * dn * PI * 0.5);
      break;
  }
}

std::complex<double> FlatVibrationNoiseSource::C_L(uint32_t n, double x) const {
  std::complex<double> CL = C_R(n, x);

  if (n == 0) {
    CL -= 2. * PI;
  }

  return CL;
}

bool FlatVibrationNoiseSource::negative_material(double x) const {
  return x < x0_ ? true : false;
}

bool FlatVibrationNoiseSource::is_inside(const Position& r) const {
  if (r.x() > low_.x() && r.y() > low_.y() && r.z() > low_.z() &&
      r.x() < hi_.x() && r.y() < hi_.y() && r.z() < hi_.z())
    return true;
  return false;
}

double FlatVibrationNoiseSource::get_x(const Position& r) const {
  switch (basis_) {
    case Basis::X:
      return r.x();
      break;

    case Basis::Y:
      return r.y();
      break;

    case Basis::Z:
      return r.z();
      break;
  }

  // Never gets here !!
  return 0.;
}

std::complex<double> FlatVibrationNoiseSource::dEt(const Position& r, double E,
                                                   double w) const {
  const double x = get_x(r);
  const bool use_neg = negative_material(x);
  const double D_Et = Delta_Et(E);

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = ((n * w0_) - w) / w;
  if (std::abs(err) > 0.01) {
    // We don't have an actual multiple of the frequency. No noise component.
    return {0., 0.};
  }

  uint32_t un = static_cast<uint32_t>(std::abs(n));

  std::complex<double> C = use_neg ? C_L(un, x) : C_R(un, x);

  if (n < 0) {
    return D_Et * C * std::exp(i * static_cast<double>(un) * PI);
  }

  return D_Et * C;
}

std::complex<double> FlatVibrationNoiseSource::dEt_Et(const Position& r,
                                                      double E,
                                                      double w) const {
  const double x = get_x(r);
  const bool use_neg = negative_material(x);
  const double D_Et = Delta_Et(E);
  const double xs = Et(x, E);
  if (xs == 0.) return {0., 0.};

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = ((n * w0_) - w) / w;
  if (std::abs(err) > 0.01) {
    // We don't have an actual multiple of the frequency. No noise component.
    return {0., 0.};
  }

  uint32_t un = static_cast<uint32_t>(std::abs(n));

  std::complex<double> C = use_neg ? C_L(un, x) : C_R(un, x);

  if (n < 0) {
    return D_Et * C * std::exp(i * static_cast<double>(un) * PI) / xs;
  }

  return D_Et * C / xs;
}

std::complex<double> FlatVibrationNoiseSource::dN(const Position& r,
                                                  uint32_t nuclide_id,
                                                  double w) const {
  const double x = get_x(r);
  const bool use_neg = negative_material(x);
  double D_N = 0.;

  // Look-up nuclide_id to get pre-computed delta
  if (Delta_N.find(nuclide_id) != Delta_N.end()) {
    D_N = Delta_N.at(nuclide_id);
  } else {
    // D_N here would be zero, so we just return.
    return {0., 0.};
  }

  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = ((n * w0_) - w) / w;
  if (std::abs(err) > 0.01) {
    // We don't have an actual multiple of the frequency. No noise component.
    return {0., 0.};
  }

  uint32_t un = static_cast<uint32_t>(std::abs(n));

  std::complex<double> C = use_neg ? C_L(un, x) : C_R(un, x);

  if (n < 0) {
    return D_N * C * std::exp(i * static_cast<double>(un) * PI);
  }

  return D_N * C;
}

std::shared_ptr<VibrationNoiseSource> make_flat_vibration_noise_source(
    const YAML::Node& snode) {
  // Get low
  if (!snode["low"] || !snode["low"].IsSequence() ||
      !(snode["low"].size() == 3)) {
    fatal_error("No valid low entry for flat vibration noise source.");
  }

  double xl = snode["low"][0].as<double>();
  double yl = snode["low"][1].as<double>();
  double zl = snode["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!snode["hi"] || !snode["hi"].IsSequence() || !(snode["hi"].size() == 3)) {
    fatal_error("No valid hi entry for flat vibration noise source.");
  }

  double xh = snode["hi"][0].as<double>();
  double yh = snode["hi"][1].as<double>();
  double zh = snode["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get frequency
  if (!snode["angular-frequency"] || !snode["angular-frequency"].IsScalar()) {
    fatal_error(
        "No valid angular-frequency entry for flat vibration noise source.");
  }

  double w0 = snode["angular-frequency"].as<double>();

  if (w0 <= 0.) {
    fatal_error(
        "Angular frequency can not be negative or zero for flat vibration "
        "noise source.");
  }

  // Get the Basis
  FlatVibrationNoiseSource::Basis basis = FlatVibrationNoiseSource::Basis::X;
  if (!snode["direction"] || !snode["direction"].IsScalar()) {
    fatal_error(
        "No valid \"direction\" entry for flat vibration noise source.");
  }
  if (snode["direction"].as<std::string>() == "x" ||
      snode["direction"].as<std::string>() == "X") {
    basis = FlatVibrationNoiseSource::Basis::X;
  } else if (snode["direction"].as<std::string>() == "y" ||
             snode["direction"].as<std::string>() == "Y") {
    basis = FlatVibrationNoiseSource::Basis::Y;
  } else if (snode["direction"].as<std::string>() == "z" ||
             snode["direction"].as<std::string>() == "Z") {
    basis = FlatVibrationNoiseSource::Basis::Z;
  } else {
    fatal_error("Invalid direction for flat vibration noise source.");
  }

  // Get the positive material
  uint32_t pos_mat_id = 0;
  if (!snode["positive-material"] || !snode["positive-material"].IsScalar()) {
    fatal_error(
        "No valid positive-material entry given for flat vibration noise "
        "source.");
  }
  pos_mat_id = snode["positive-material"].as<uint32_t>();

  std::shared_ptr<Material> pos_mat = nullptr;
  if (materials.find(pos_mat_id) == materials.end()) {
    std::stringstream mssg;
    mssg << "poitive-material with id " << pos_mat_id
         << " not found for flat vibration noise source.";
    fatal_error(mssg.str());
  }
  pos_mat = materials[pos_mat_id];

  // Get the negative material
  uint32_t neg_mat_id = 0;
  if (!snode["negative-material"] || !snode["negative-material"].IsScalar()) {
    fatal_error(
        "No valid negative-material entry given for flat vibration noise "
        "source.");
  }
  neg_mat_id = snode["negative-material"].as<uint32_t>();

  std::shared_ptr<Material> neg_mat = nullptr;
  if (materials.find(neg_mat_id) == materials.end()) {
    std::stringstream mssg;
    mssg << "negative-material with id " << pos_mat_id
         << " not found for flat vibration noise source.";
    fatal_error(mssg.str());
  }
  neg_mat = materials[neg_mat_id];

  return std::make_shared<FlatVibrationNoiseSource>(r_low, r_hi, basis, pos_mat,
                                                    neg_mat, w0);
}
