#include <materials/material_helper.hpp>
#include <noise_source/propagating_channel_oscillation_noise_source.hpp>
#include <noise_source/square_oscillation_noise_source.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/direction.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <sstream>

PropagatingChannelOscillationNoiseSource::
    PropagatingChannelOscillationNoiseSource(Position low, Position hi,
                                             double radius, char axis,
                                             double eps_tot, double eps_fis,
                                             double eps_sct,
                                             double angular_frequency,
                                             double velocity, double phase)
    : low_(low),
      hi_(hi),
      r_origin_(0.5 * (low_.x() + hi_.x()), 0.5 * (low_.y() + hi_.y()),
                0.5 * (low_.z() + hi_.z())),
      radius_(radius),
      axis_(axis),
      w0_(angular_frequency),
      phase_(phase),
      velocity_(velocity),
      eps_t_(eps_tot),
      eps_f_(eps_fis),
      eps_s_(eps_sct) {
  // Check low and high
  if (low_.x() >= hi_.x() || low_.y() >= hi_.y() || low_.z() >= hi_.z()) {
    fatal_error("Low is greater than or equal to hi.");
  }

  // Check radius
  if (radius_ <= 0.) {
    fatal_error("Radius is <= 0.");
  }

  // Make sure axis is valid
  if (axis_ != 'x' && axis_ != 'y' && axis_ != 'z') {
    std::stringstream mssg;
    mssg << "Unknown axis \"" << axis_ << "\".";
    fatal_error(mssg.str());
  }

  // Check velocity
  if (velocity_ == 0.) {
    fatal_error("Negative propagation velocity.");
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

bool PropagatingChannelOscillationNoiseSource::is_inside(
    const Position& r) const {
  if (r.x() > low_.x() && r.y() > low_.y() && r.z() > low_.z() &&
      r.x() < hi_.x() && r.y() < hi_.y() && r.z() < hi_.z()) {
    // Now we check the radius
    const double rad = calc_center_line_radius(r);

    if (rad > radius_) return true;

    return false;
  }

  return false;
}

std::complex<double> PropagatingChannelOscillationNoiseSource::dEt(
    const Position& r, double E, double w) const {
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

  const double tot_phase = calc_total_phase(r);

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * xs * PI) *
           std::exp(std::complex<double>(0., tot_phase));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * xs * PI) *
           std::exp(std::complex<double>(0., -tot_phase));
  } else {
    return {0., 0.};
  }
}

std::complex<double> PropagatingChannelOscillationNoiseSource::dEt_Et(
    const Position& r, double /*E*/, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  const double tot_phase = calc_total_phase(r);

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * PI) *
           std::exp(std::complex<double>(0., tot_phase));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_t_ * PI) *
           std::exp(std::complex<double>(0., -tot_phase));
  } else {
    return {0., 0.};
  }
}

std::complex<double> PropagatingChannelOscillationNoiseSource::dEf_Ef(
    const Position& r, double /*E*/, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  const double tot_phase = calc_total_phase(r);

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_ * PI) *
           std::exp(std::complex<double>(0., tot_phase));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_f_ * PI) *
           std::exp(std::complex<double>(0., -tot_phase));
  } else {
    return {0., 0.};
  }
}

std::complex<double>
PropagatingChannelOscillationNoiseSource::dEelastic_Eelastic(const Position& r,
                                                             double /*E*/,
                                                             double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  const double tot_phase = calc_total_phase(r);

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., tot_phase));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., -tot_phase));
  } else {
    return {0., 0.};
  }
}

std::complex<double> PropagatingChannelOscillationNoiseSource::dEmt_Emt(
    uint32_t /*mt*/, const Position& r, double /*E*/, double w) const {
  // Get the frequency multiple n
  int32_t n = static_cast<int32_t>(std::round(w / w0_));

  double err = (n * w0_ - w) / w;

  const double tot_phase = calc_total_phase(r);

  if (n == 1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., tot_phase));
  } else if (n == -1 && std::abs(err) < 0.01) {
    return std::complex<double>(0., -eps_s_ * PI) *
           std::exp(std::complex<double>(0., -tot_phase));
  } else {
    return {0., 0.};
  }
}

double PropagatingChannelOscillationNoiseSource::calc_center_line_radius(
    const Position& r) const {
  if (axis_ == 'x') {
    const double dy = r.y() - r_origin_.y();
    const double dz = r.z() - r_origin_.z();
    return std::sqrt(dy * dy + dz * dz);
  } else if (axis_ == 'y') {
    const double dx = r.x() - r_origin_.x();
    const double dz = r.z() - r_origin_.z();
    return std::sqrt(dx * dx + dz * dz);
  } else {
    const double dx = r.x() - r_origin_.x();
    const double dy = r.y() - r_origin_.y();
    return std::sqrt(dy * dy + dx * dx);
  }
}

double PropagatingChannelOscillationNoiseSource::calc_total_phase(
    const Position& r) const {
  if (axis_ == 'x') {
    return phase_ - (w0_ * r.x() / velocity_);
  } else if (axis_ == 'y') {
    return phase_ - (w0_ * r.y() / velocity_);
  } else {
    return phase_ - (w0_ * r.z() / velocity_);
  }
}

std::shared_ptr<OscillationNoiseSource>
make_propagating_channel_oscillation_noise_source(const YAML::Node& snode) {
  // Get low
  if (!snode["low"] || !snode["low"].IsSequence() ||
      !(snode["low"].size() == 3)) {
    fatal_error(
        "No valid low entry for propagating channel oscillation noise source.");
  }

  double xl = snode["low"][0].as<double>();
  double yl = snode["low"][1].as<double>();
  double zl = snode["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!snode["hi"] || !snode["hi"].IsSequence() || !(snode["hi"].size() == 3)) {
    fatal_error(
        "No valid hi entry for propagating channel oscillation noise source.");
  }

  double xh = snode["hi"][0].as<double>();
  double yh = snode["hi"][1].as<double>();
  double zh = snode["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get radius
  if (!snode["radius"] || !snode["radius"].IsScalar()) {
    fatal_error(
        "No valid radius entry for propagating channel oscillation noise "
        "source.");
  }
  const double rad = snode["radius"].as<double>();
  if (rad <= 0.) {
    fatal_error(
        "Noise source of type propagating channel oscillation has a radius "
        "that is <= 0.");
  }

  // Get axis
  if (!snode["axis"] || !snode["axis"].IsScalar()) {
    fatal_error(
        "No valid axis entry for propagating channel oscillation noise "
        "source.");
  }
  std::string axis_str = snode["axis"].as<std::string>();
  if (axis_str != "x" && axis_str != "y" && axis_str != "z") {
    std::stringstream mssg;
    mssg << "Noise source of type propagating channel oscillation was provided "
            "unkown axis \""
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
        "No valid angular-frequency entry for propagating channel oscillation "
        "noise source.");
  }

  double w0 = snode["angular-frequency"].as<double>();

  if (w0 <= 0.) {
    fatal_error(
        "Angular frequency can not be negative or zero for propagating channel "
        "oscillation noise source.");
  }

  // Get phase (is present)
  double phase = 0.;
  if (snode["phase"] && snode["phase"].IsScalar()) {
    phase = snode["phase"].as<double>();
  } else if (snode["phase"]) {
    fatal_error(
        "Invalid phase entry for propagating channel oscillation noise "
        "source.");
  }

  // Get velocity
  double velocity = 0.;
  if (!snode["velocity"] || snode["velocity"].IsScalar() == false) {
    fatal_error(
        "No valid velocity entry in propagating channel oscillation noise "
        "source.");
  }
  velocity = snode["velocity"].as<double>();
  if (velocity == 0.) {
    fatal_error(
        "Provided velocity in propagating channel oscillation noise source is "
        "zero.");
  }

  // Get epsilons
  if (!snode["epsilon-total"] || !snode["epsilon-total"].IsScalar()) {
    fatal_error(
        "No valid epsilon-total entry for propagating channel oscillation "
        "noise source.");
  }

  double eps_t = snode["epsilon-total"].as<double>();

  if (eps_t <= 0.) {
    fatal_error(
        "Epsilon total can not be negative or zero for propagating channel "
        "oscillation noise source.");
  }

  if (!snode["epsilon-fission"] || !snode["epsilon-fission"].IsScalar()) {
    fatal_error(
        "No valid epsilon-fission entry for propagating channel oscillation "
        "noise source.");
  }

  double eps_f = snode["epsilon-fission"].as<double>();

  if (eps_f <= 0.) {
    fatal_error(
        "Epsilon fission can not be negative or zero for propagating channel "
        "oscillation noise source.");
  }

  if (!snode["epsilon-scatter"] || !snode["epsilon-scatter"].IsScalar()) {
    fatal_error(
        "No valid epsilon-scatter entry for propagating channel oscillation "
        "noise source.");
  }

  double eps_s = snode["epsilon-scatter"].as<double>();

  if (eps_s <= 0.) {
    fatal_error(
        "Epsilon scatter can not be negative or zero for propagating channel "
        "oscillation noise source.");
  }

  return std::make_shared<PropagatingChannelOscillationNoiseSource>(
      r_low, r_hi, rad, axis, eps_t, eps_f, eps_s, w0, velocity, phase);
}
