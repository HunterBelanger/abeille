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
#include <noise_source/cylindrical_oscillation_noise_source.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/direction.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

#include <sstream>

CylindricalOscillationNoiseSource::CylindricalOscillationNoiseSource(
    Position origin, double len, double rad, char ax,
    const pndl::Tabulated1D& eps_tot, const pndl::Tabulated1D& eps_fis,
    const pndl::Tabulated1D& eps_sct, double angular_frequency, double phase)
    : origin_(origin),
      length_(len),
      radius_(rad),
      axis_(ax),
      w0_(angular_frequency),
      phase_(phase),
      eps_t_(eps_tot),
      eps_f_(eps_fis),
      eps_s_(eps_sct) {
  // Check length and radius
  if (length_ <= 0.) {
    fatal_error("Length is <= 0.");
  }

  if (radius_ <= 0.) {
    fatal_error("Radius is <= 0.");
  }

  // Make sure axis is valid
  if (axis_ != 'x' && axis_ != 'y' && axis_ != 'z') {
    std::stringstream mssg;
    mssg << "Unknown axis \"" << axis_ << "\".";
    fatal_error(mssg.str());
  }

  // Make sure frequency is positive
  if (w0_ <= 0.) {
    fatal_error("Negative or zero noise frequency.");
  }

  // Make sure epsilons are positive
  for (const auto et : eps_t_.y()) {
    if (et < 0.) {
      fatal_error("Negative epsilon total.");
    }
  }

  for (const auto ef : eps_f_.y()) {
    if (ef < 0.) {
      fatal_error("Negative epsilon fission.");
    }
  }

  for (const auto es : eps_s_.y()) {
    if (es < 0.) {
      fatal_error("Negative epsilon scatter.");
    }
  }
}

bool CylindricalOscillationNoiseSource::is_inside(const Position& r) const {
  const Position r_loc = r - origin_;

  if (axis_ == 'x') {
    if (std::abs(r_loc.x()) >= 0.5 * length_) {
      return false;
    }

    if (std::abs(r_loc.y()) >= radius_ || std::abs(r_loc.z()) >= radius_) {
      return false;
    }

    if (std::sqrt(r_loc.y() * r_loc.y() + r_loc.z() * r_loc.z()) >= radius_) {
      return false;
    }

    return true;
  } else if (axis_ == 'y') {
    if (std::abs(r_loc.y()) >= 0.5 * length_) {
      return false;
    }

    if (std::abs(r_loc.x()) >= radius_ || std::abs(r_loc.z()) >= radius_) {
      return false;
    }

    if (std::sqrt(r_loc.x() * r_loc.x() + r_loc.z() * r_loc.z()) >= radius_) {
      return false;
    }

    return true;
  } else {
    if (std::abs(r_loc.z()) >= 0.5 * length_) {
      return false;
    }

    if (std::abs(r_loc.y()) >= radius_ || std::abs(r_loc.x()) >= radius_) {
      return false;
    }

    if (std::sqrt(r_loc.y() * r_loc.y() + r_loc.x() * r_loc.x()) >= radius_) {
      return false;
    }

    return true;
  }
}

std::complex<double> CylindricalOscillationNoiseSource::dEt(const Position& r,
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
    return std::complex<double>(0., -eps_t_(E) * xs * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_(E) * xs * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> CylindricalOscillationNoiseSource::dEt_Et(
    const Position& /*r*/, double E, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_(E) * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_(E) * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> CylindricalOscillationNoiseSource::dEf_Ef(
    const Position& /*r*/, double E, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_(E) * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_(E) * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> CylindricalOscillationNoiseSource::dEelastic_Eelastic(
    const Position& /*r*/, double E, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_(E) * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_(E) * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::complex<double> CylindricalOscillationNoiseSource::dEmt_Emt(
    uint32_t /*mt*/, const Position& /*r*/, double E, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_(E) * PI) *
           std::exp(std::complex<double>(0., phase_));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_(E) * PI) *
           std::exp(std::complex<double>(0., -phase_));
  } else {
    return {0., 0.};
  }
}

std::shared_ptr<OscillationNoiseSource>
make_cylindrical_oscillation_noise_source(const YAML::Node& snode) {
  // Get origin
  if (!snode["origin"] || !snode["origin"].IsSequence() ||
      !(snode["origin"].size() == 3)) {
    fatal_error(
        "No valid origin entry for cylindrical oscillation noise source.");
  }

  const double xo = snode["origin"][0].as<double>();
  const double yo = snode["origin"][1].as<double>();
  const double zo = snode["origin"][2].as<double>();
  const Position r_orig(xo, yo, zo);

  // Get length
  if (!snode["length"] || !snode["length"].IsScalar()) {
    fatal_error(
        "No valid length entry for cylindrical oscillation noise source.");
  }
  const double len = snode["length"].as<double>();
  if (len <= 0.) {
    fatal_error(
        "Noise source of type cylindrical oscillation has a length that is <= "
        "0.");
  }

  // Get radius
  if (!snode["radius"] || !snode["radius"].IsScalar()) {
    fatal_error(
        "No valid radius entry for cylindrical oscillation noise source.");
  }
  const double rad = snode["radius"].as<double>();
  if (rad <= 0.) {
    fatal_error(
        "Noise source of type cylindrical oscillation has a radius that is <= "
        "0.");
  }

  // Get axis
  if (!snode["axis"] || !snode["axis"].IsScalar()) {
    fatal_error(
        "No valid axis entry for cylindrical oscillation noise source.");
  }
  std::string axis_str = snode["axis"].as<std::string>();
  if (axis_str != "x" && axis_str != "y" && axis_str != "z") {
    std::stringstream mssg;
    mssg << "Noise source of type cylindrical oscillation was provided unkown "
            "axis \""
         << axis_str << "\".";
    fatal_error(mssg.str());
  }
  char axis;
  if (axis_str == "x")
    axis = 'x';
  else if (axis_str == "y")
    axis = 'y';
  else
    axis = 'z';

  // Get frequency
  if (!snode["angular-frequency"] || !snode["angular-frequency"].IsScalar()) {
    fatal_error(
        "No valid angular-frequency entry for cylindrical oscillation noise "
        "source.");
  }

  double w0 = snode["angular-frequency"].as<double>();

  if (w0 <= 0.) {
    fatal_error(
        "Angular frequency can not be negative or zero for oscillation noise "
        "source.");
  }

  // Get phase (is present)
  double phase = 0.;
  if (snode["phase"] && snode["phase"].IsScalar()) {
    phase = snode["phase"].as<double>();
  } else if (snode["phase"]) {
    fatal_error(
        "Invalid phase entry for cylindrical oscillation noise source.");
  }

  // Get epsilons
  if (!snode["epsilon-total"] || !snode["epsilon-total"].IsMap()) {
    fatal_error("No valid epsilon-total entry for oscillation noise source.");
  }
  pndl::Tabulated1D eps_t =
      make_tabulated_1d(snode["epsilon-total"], "energy", "value");
  for (const auto et : eps_t.y()) {
    if (et < 0.) {
      fatal_error(
          "Epsilon total can not be negative for oscillation noise source.");
    }
  }
  if (eps_t.min_x() > settings::min_energy ||
      eps_t.max_x() < settings::max_energy) {
    fatal_error(
        "Epsilon total energies are smaller than min and max particle "
        "energies.");
  }

  if (!snode["epsilon-fission"] || !snode["epsilon-fission"].IsMap()) {
    fatal_error("No valid epsilon-fission entry for oscillation noise source.");
  }
  pndl::Tabulated1D eps_f =
      make_tabulated_1d(snode["epsilon-fission"], "energy", "value");
  for (const auto ef : eps_f.y()) {
    if (ef < 0.) {
      fatal_error(
          "Epsilon fission can not be negative for oscillation noise source.");
    }
  }
  if (eps_f.min_x() > settings::min_energy ||
      eps_f.max_x() < settings::max_energy) {
    fatal_error(
        "Epsilon fission energies are smaller than min and max particle "
        "energies.");
  }

  if (!snode["epsilon-scatter"] || !snode["epsilon-scatter"].IsMap()) {
    fatal_error("No valid epsilon-scatter entry for oscillation noise source.");
  }
  pndl::Tabulated1D eps_s =
      make_tabulated_1d(snode["epsilon-scatter"], "energy", "value");
  for (const auto es : eps_s.y()) {
    if (es < 0.) {
      fatal_error(
          "Epsilon scatter can not be negative for oscillation noise source.");
    }
  }
  if (eps_s.min_x() > settings::min_energy ||
      eps_s.max_x() < settings::max_energy) {
    fatal_error(
        "Epsilon scatter energies are smaller than min and max particle "
        "energies.");
  }

  return std::make_shared<CylindricalOscillationNoiseSource>(
      r_orig, len, rad, axis, eps_t, eps_f, eps_s, w0, phase);
}
