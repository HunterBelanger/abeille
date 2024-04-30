#include <tallies/box_position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/mesh_position_filter.hpp>
#include <utils/constants.hpp>

bool CartesianFilter::find_entry_point(Position& r, const Direction& u,
                                       double& d_flight) const {
  const double ux_inv = 1. / u.x();
  const double uy_inv = 1. / u.y();
  const double uz_inv = 1. / u.z();

  double d_min = (r_low.x() - r.x()) * ux_inv;
  double d_max = (r_high.x() - r.x()) * ux_inv;

  if (d_min > d_max) {
    std::swap(d_min, d_max);
  }

  double d_y_min = (r_low.y() - r.y()) * uy_inv;
  double d_y_max = (r_high.y() - r.y()) * uy_inv;

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

  double d_z_min = (r_low.z() - r.z()) * uz_inv;
  double d_z_max = (r_high.z() - r.z()) * uz_inv;

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

void CartesianFilter::initialize_indices(const Position& r, const Direction& u,
                                         int& i, int& j, int& k,
                                         std::array<int, 3>& on) {
  i = static_cast<int>(std::floor((r.x() - r_low.x()) * dx_inv));
  j = static_cast<int>(std::floor((r.y() - r_low.y()) * dy_inv));
  k = static_cast<int>(std::floor((r.z() - r_low.z()) * dz_inv));
  on.fill(0);

  // Must handle case of being on a tile boundary
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get tile boundaries
  double xl = xc - 0.5 * dx;
  double xh = xc + 0.5 * dx;
  double yl = yc - 0.5 * dy;
  double yh = yc + 0.5 * dy;
  double zl = zc - 0.5 * dz;
  double zh = zc + 0.5 * dz;

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

void CartesianFilter::update_indices(int key, int& i, int& j, int& k,
                                     std::array<int, 3>& on) {
  // Must initially fill with zero, so that we don't stay on top
  // of other surfaces the entire time
  on.fill(0);

  switch (key) {
    case -1:
      i--;
      on[0] = 1;
      break;

    case 1:
      i++;
      on[0] = -1;
      break;

    case -2:
      j--;
      on[1] = 1;
      break;

    case 2:
      j++;
      on[1] = -1;
      break;

    case -3:
      k--;
      on[2] = 1;
      break;

    case 3:
      k++;
      on[2] = -1;
      break;

    default:
      break;
  }
}

std::pair<double, int> CartesianFilter::distance_to_next_index(
    const Position& r, const Direction& u, const std::array<int, 3>& on, int i,
    int j, int k) {
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get relative position in cell
  Position r_tile(r.x() - xc, r.y() - yc, r.z() - zc);

  // Set our initial value for the distance and the index change
  double dist = INF;
  int key = 0;

  // Check all six sides
  const double diff_xl = -dx * 0.5 - r_tile.x();
  const double diff_xh = dx * 0.5 - r_tile.x();
  const double diff_yl = -dy * 0.5 - r_tile.y();
  const double diff_yh = dy * 0.5 - r_tile.y();
  const double diff_zl = -dz * 0.5 - r_tile.z();
  const double diff_zh = dz * 0.5 - r_tile.z();

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

// make the cartesian filter
std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node) {
  std::shared_ptr<CartesianFilter> cartesian_filter_ = nullptr;

  if (!node["Position-Filter"]) {
    fatal_error("Position-Filter is not given.");
  }
  const std::string cartesian_filter_type =
      node["Position-Filter"].as<std::string>();

  if (cartesian_filter_type == "box") {
    cartesian_filter_ = make_box_position_filter<CartesianFilter>(node);
  }

  if (cartesian_filter_type == "regular-reactangular-mesh") {
    cartesian_filter_ = make_mesh_position_filter<CartesianFilter>(node);
  }

  return cartesian_filter_;
}