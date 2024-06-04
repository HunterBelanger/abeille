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
    std::stringstream mssg;
    mssg << "Negative propagation velocity in oscillation noise source.";
    fatal_error(mssg.str());
  }

  // Make sure frequency is positive
  if (w0_ <= 0.) {
    fatal_error("Negative or zero frequency oscillation frequency.");
  }

  // Make sure epsilons are positive
  if (eps_t_ <= 0.) {
    fatal_error("Negative or zero epsilon total in oscillation noise source.");
  }

  if (eps_f_ <= 0.) {
    fatal_error(
        "Negative or zero epsilon fission in oscillation noise source.");
  }

  if (eps_s_ <= 0.) {
    fatal_error(
        "Negative or zero epsilon scatter in oscillation noise source.");
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