#include <tallies/box_position_filter.hpp>
#include <utils/constants.hpp>
#include <utils/output.hpp>

BoxPositionFilter::BoxPositionFilter(const Position r_low,
                                     const Position r_high)
    : CartesianFilter(r_low, r_high),
      dx_(),
      dy_(),
      dz_(),
      dx_inv_(),
      dy_inv_(),
      dz_inv_() {
  if ((r_low_.x() > r_high_.x()) || (r_low_.y() > r_high_.y()) ||
      (r_low_.z() > r_high_.z()))
    fatal_error(
        " Corrdinates of \"low\" position are higher than \"high\" "
        "position.\n");

  dx_ = (r_high_.x() - r_low_.x());
  dy_ = (r_high_.y() - r_low_.y());
  dz_ = (r_high_.z() - r_low_.z());

  dx_inv_ = 1. / dx_;
  dy_inv_ = 1. / dy_;
  dz_inv_ = 1. / dz_;
}

std::vector<TracklengthDistance> BoxPositionFilter::get_indices_tracklength(
    const Tracker& trkr, double d_flight) {
  std::vector<TracklengthDistance> indexes_tracklength;
  TracklengthDistance trlen_d;

  Position r = trkr.r();
  const Direction u_ = trkr.u();
  bool inside_bin = false;

  int i = 0, j = 0, k = 0;
  std::array<int, 3> on;
  on.fill(0.0);

  // Check if particle inside the bin
  if ((r_low_.x() <= r.x() && r_high_.x() >= r.x()) &&
      (r_low_.y() <= r.y() && r_high_.y() >= r.y()) &&
      (r_low_.z() <= r.z() && r_high_.z() >= r.z())) {
    inside_bin = true;
  }

  // it is possible that particle is not inside the box, but can intersect
  if (inside_bin == false) {
    if (find_entry_point(r, u_, d_flight) == false) {
      return indexes_tracklength;
    }
    initialize_indices(r, u_, i, j, k, on);
    if (i == 0 && j == 0 && k == 0) {
      inside_bin = true;
    } else {
      // This is a problem, in theory, we should now be inside the tally
      // region. We will therefore spew a warning here.
      warning("Could not locate tile after fast forward to mesh entry.\n");
    }
  }

  // Distance remaining to tally
  double distance_remaining = d_flight;

  auto next_tile = distance_to_next_index(r, u_, on, i, j, k);

  if (next_tile.first == INF) {
    // Something went wrong.... Don't score.
    Output::instance().save_warning(
        "Problem encountered with mesh tally in box_tally.\n");
    return indexes_tracklength;

  } else if (next_tile.first < 0.) {
    // Something went wrong.... Don't score.
    warning("Negative distance encountered with mesh tally.");
  }

  double d_tile = std::min(next_tile.first, distance_remaining);

  // Make the score if we are in a valid cell
  if (i == 0 && j == 0 && k == 0) {
    trlen_d.indexes_ = StaticVector3{0};
    trlen_d.distance_in_bin = d_tile;
    indexes_tracklength.push_back(trlen_d);
  }
  return indexes_tracklength;
}

bool BoxPositionFilter::find_entry_point(Position& r, const Direction& u,
                                         double& d_flight) const {
  const double ux_inv = 1. / u.x();
  const double uy_inv = 1. / u.y();
  const double uz_inv = 1. / u.z();

  double d_min = (r_low_.x() - r.x()) * ux_inv;
  double d_max = (r_high_.x() - r.x()) * ux_inv;

  if (d_min > d_max) {
    std::swap(d_min, d_max);
  }

  double d_y_min = (r_low_.y() - r.y()) * uy_inv;
  double d_y_max = (r_high_.y() - r.y()) * uy_inv;

  if (d_y_min > d_y_max) {
    std::swap(d_y_min, d_y_max);
  }

  if ((d_min > d_y_max) || (d_y_min > d_max)) {
    return false;
  }

  if (d_y_min > d_min) {
    d_min = d_y_min;
  }

  if (d_y_max < d_max) {
    d_max = d_y_max;
  }

  double d_z_min = (r_low_.z() - r.z()) * uz_inv;
  double d_z_max = (r_high_.z() - r.z()) * uz_inv;

  if (d_z_min > d_z_max) {
    std::swap(d_z_min, d_z_max);
  }

  if ((d_min > d_z_max) || (d_z_min > d_max)) {
    return false;
  }

  if (d_z_min > d_min) {
    d_min = d_z_min;
  }

  if (d_z_max < d_max) {
    d_max = d_z_max;
  }

  if (d_max < d_min) {
    std::swap(d_max, d_min);
  }

  if ((d_max < 0.) && (d_min < 0.)) {
    return false;
  }

  if (d_min < 0.) {
    // If we are here, this means that r is actually inside the mesh, but is
    // really close to the edge, and we have a direction taking us out.
    // We should return false here, so that we don't score anything for this
    // particle track.
    return false;
  }

  // If we get here, we intersect the box. Let's update the position and the
  // flight distance.
  r = r + d_min * u;
  d_flight -= d_min;
  return true;
}

void BoxPositionFilter::initialize_indices(const Position& r,
                                           const Direction& u, int& i, int& j,
                                           int& k, std::array<int, 3>& on) {
  i = static_cast<int>(std::floor((r.x() - r_low_.x()) * dx_inv_));
  j = static_cast<int>(std::floor((r.y() - r_low_.y()) * dy_inv_));
  k = static_cast<int>(std::floor((r.z() - r_low_.z()) * dz_inv_));
  on.fill(0);

  // Must handle case of being on a tile boundary
  // Get position at center of current tile
  double xc = r_low_.x() + i * dx_ + 0.5 * dx_;
  double yc = r_low_.y() + j * dy_ + 0.5 * dy_;
  double zc = r_low_.z() + k * dz_ + 0.5 * dz_;

  // Get tile boundaries
  double xl = xc - 0.5 * dx_;
  double xh = xc + 0.5 * dx_;
  double yl = yc - 0.5 * dy_;
  double yh = yc + 0.5 * dy_;
  double zl = zc - 0.5 * dz_;
  double zh = zc + 0.5 * dz_;

  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      i--;
      on[0] = 1;
    } else {
      on[0] = -1;
    }
  } else if (std::abs(xh - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      on[0] = 1;
    } else {
      i++;
      on[0] = -1;
    }
  }

  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      j--;
      on[1] = 1;
    } else {
      on[1] = -1;
    }
  } else if (std::abs(yh - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      on[1] = 1;
    } else {
      j++;
      on[1] = -1;
    }
  }

  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      k--;
      on[2] = 1;
    } else {
      on[2] = -1;
    }
  } else if (std::abs(zh - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      on[2] = 1;
    } else {
      k++;
      on[2] = -1;
    }
  }
}

std::pair<double, int> BoxPositionFilter::distance_to_next_index(
    const Position& r, const Direction& u, const std::array<int, 3>& on, int i,
    int j, int k) {
  // Get position at center of current tile
  double xc = r_low_.x() + i * dx_ + 0.5 * dx_;
  double yc = r_low_.y() + j * dy_ + 0.5 * dy_;
  double zc = r_low_.z() + k * dz_ + 0.5 * dz_;

  // Get relative position in cell
  Position r_tile(r.x() - xc, r.y() - yc, r.z() - zc);

  // Set our initial value for the distance and the index change
  double dist = INF;
  int key = 0;

  // Check all six sides
  const double diff_xl = -dx_ * 0.5 - r_tile.x();
  const double diff_xh = dx_ * 0.5 - r_tile.x();
  const double diff_yl = -dy_ * 0.5 - r_tile.y();
  const double diff_yh = dy_ * 0.5 - r_tile.y();
  const double diff_zl = -dz_ * 0.5 - r_tile.z();
  const double diff_zh = dz_ * 0.5 - r_tile.z();

  const double ux_inv = 1. / u.x();
  const double uy_inv = 1. / u.y();
  const double uz_inv = 1. / u.z();

  const double d_xl = diff_xl * ux_inv;
  const double d_xh = diff_xh * ux_inv;
  const double d_yl = diff_yl * uy_inv;
  const double d_yh = diff_yh * uy_inv;
  const double d_zl = diff_zl * uz_inv;
  const double d_zh = diff_zh * uz_inv;

  if (d_xl > 0. && d_xl < dist && on[0] != -1) {
    dist = d_xl;
    key = -1;
  }

  if (d_xh > 0. && d_xh < dist && on[0] != 1) {
    dist = d_xh;
    key = 1;
  }

  if (d_yl > 0. && d_yl < dist && on[1] != -1) {
    dist = d_yl;
    key = -2;
  }

  if (d_yh > 0. && d_yh < dist && on[1] != 1) {
    dist = d_yh;
    key = 2;
  }

  if (d_zl > 0. && d_zl < dist && on[2] != -1) {
    dist = d_zl;
    key = -3;
  }

  if (d_zh > 0. && d_zh < dist && on[2] != 1) {
    dist = d_zh;
    key = 3;
  }

  return {dist, key};
}

// make the cartesian filter or position filter
std::shared_ptr<BoxPositionFilter> make_box_position_filter(
    const YAML::Node& node) {
  if (!node["low"])
    fatal_error(
        "For box position-filter \"low\" co-ordinates is not provided.");
  if (!node["high"])
    fatal_error(
        "For box position-filter \"high\" co-ordinates is not provided.");

  std::vector<double> low_point = node["low"].as<std::vector<double>>();
  std::vector<double> high_point = node["high"].as<std::vector<double>>();

  Position r_low(low_point[0], low_point[1], low_point[2]);
  Position r_high(high_point[0], high_point[1], high_point[2]);

  std::shared_ptr<BoxPositionFilter> box_type_filter =
      std::make_shared<BoxPositionFilter>(r_low, r_high);

  return box_type_filter;
}