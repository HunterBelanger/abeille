#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/constants.hpp>
#include <utils/output.hpp>

#include <sstream>

RegularCartesianMeshFilter::RegularCartesianMeshFilter(Position r_low,
                                                       Position r_high,
                                                       size_t nx_, size_t ny_,
                                                       size_t nz_,
                                                       std::size_t id)
    : CartesianFilter(r_low, r_high, id),
      dx_(),
      dy_(),
      dz_(),
      dx_inv_(),
      dy_inv_(),
      dz_inv_(),
      Nx_(nx_),
      Ny_(ny_),
      Nz_(nz_),
      x_index_(),
      y_index_(),
      z_index_() {
  if ((r_low_.x() >= r_high_.x()) || (r_low_.y() >= r_high_.y()) ||
      (r_low_.z() >= r_high_.z()))
    fatal_error("In position-filter with id: " + std::to_string(id) +
                ", coordinates of \"low\" position are >= than \"high\" "
                "position.");

  if (Nx_ == 0 || Ny_ == 0 || Nz_ == 0)
    fatal_error("In position-filter with id: " + std::to_string(id) +
                ", the number of bins in any direction must be a non-zero.");

  dx_ = (r_high_.x() - r_low_.x()) / static_cast<double>(Nx_);
  dy_ = (r_high_.y() - r_low_.y()) / static_cast<double>(Ny_);
  dz_ = (r_high_.z() - r_low_.z()) / static_cast<double>(Nz_);

  dx_inv_ = 1. / dx_;
  dy_inv_ = 1. / dy_;
  dz_inv_ = 1. / dz_;

  x_index_ = 0;
  y_index_ = 1;
  z_index_ = 2;
  if (Nx_ == 1) {
    x_index_ = 0;
    y_index_--;
    z_index_--;
  }

  if (Ny_ == 1) {
    y_index_ = 0;
    z_index_--;
  }

  if (Nz_ == 1) {
    z_index_ = 0;
  }
}

StaticVector3 RegularCartesianMeshFilter::get_indices(
    const Tracker& tktr) const {
  StaticVector3 indices;
  const Position r = tktr.r();
  const int index_x =
      static_cast<int>(std::floor((r.x() - r_low_.x()) * dx_inv_));
  const int index_y =
      static_cast<int>(std::floor((r.y() - r_low_.y()) * dy_inv_));
  const int index_z =
      static_cast<int>(std::floor((r.z() - r_low_.z()) * dz_inv_));

  if ((index_x >= 0 && index_x < static_cast<int>(Nx_)) &&
      (index_y >= 0 && index_y < static_cast<int>(Ny_)) &&
      (index_z >= 0 && index_z < static_cast<int>(Nz_))) {
    indices = reduce_dimension(index_x, index_y, index_z);
  }

  return indices;
}

StaticVector3 RegularCartesianMeshFilter::get_position_index(
    const Position& r) const {
  StaticVector3 indices;
  const int index_x =
      static_cast<int>(std::floor((r.x() - r_low_.x()) * dx_inv_));
  const int index_y =
      static_cast<int>(std::floor((r.y() - r_low_.y()) * dy_inv_));
  const int index_z =
      static_cast<int>(std::floor((r.z() - r_low_.z()) * dz_inv_));

  if ((index_x >= 0 && index_x < static_cast<int>(Nx_)) &&
      (index_y >= 0 && index_y < static_cast<int>(Ny_)) &&
      (index_z >= 0 && index_z < static_cast<int>(Nz_))) {
    indices = reduce_dimension(index_x, index_y, index_z);
  }

  return indices;
}

double RegularCartesianMeshFilter::x_min(const StaticVector3& index) const {
  if (Nx_ == 1) {
    return r_low_.x();
  }
  return r_low_.x() + static_cast<double>(index[x_index_]) * dx_;
}

double RegularCartesianMeshFilter::x_max(const StaticVector3& index) const {
  if (Nx_ == 1) {
    return r_high_.x();
  }
  return r_low_.x() + static_cast<double>(index[x_index_]) * dx_ + dx_;
}

double RegularCartesianMeshFilter::y_min(const StaticVector3& index) const {
  if (Ny_ == 1) {
    return r_low_.y();
  }

  return r_low_.y() + static_cast<double>(index[y_index_]) * dy_;
}

double RegularCartesianMeshFilter::y_max(const StaticVector3& index) const {
  if (Ny_ == 1) return r_high_.y();

  return r_low_.y() + static_cast<double>(index[y_index_]) * dy_ + dy_;
}

double RegularCartesianMeshFilter::z_min(const StaticVector3& index) const {
  if (Nz_ == 1) return r_low_.z();

  return r_low_.z() + static_cast<double>(index[z_index_]) * dz_;
}

double RegularCartesianMeshFilter::z_max(const StaticVector3& index) const {
  if (Nz_ == 1) return r_high_.z();

  return r_low_.z() + static_cast<double>(index[z_index_]) * dz_ + dz_;
}

std::vector<TracklengthDistance>
RegularCartesianMeshFilter::get_indices_tracklength(const Tracker& trkr,
                                                    double d_flight) const {
  std::vector<TracklengthDistance> indices_tracklength;
  TracklengthDistance trlen_d;

  Position r = trkr.r();
  const Direction u_ = trkr.u();
  bool inside_bin = false;

  int i = 0, j = 0, k = 0;
  std::array<int, 3> on;
  on.fill(0);
  initialize_indices(r, u_, i, j, k, on);

  // Check if particle inside the bin
  if (i >= 0 && i < static_cast<int>(Nx_) && j >= 0 &&
      j < static_cast<int>(Ny_) && k >= 0 && k < static_cast<int>(Nz_)) {
    inside_bin = true;
  }

  // it is possible that particle is not inside the box, but can intersect
  if (inside_bin == false) {
    if (find_entry_point(r, u_, d_flight) == false) {
      return indices_tracklength;
    }
    initialize_indices(r, u_, i, j, k, on);
    if (i >= 0 && i < static_cast<int>(Nx_) && j >= 0 &&
        j < static_cast<int>(Ny_) && k >= 0 && k < static_cast<int>(Nz_)) {
      inside_bin = true;
    } else {
      // This is a problem, in theory, we should now be inside the tally
      // region. We will therefore spew a warning here.
      warning("Could not locate tile after fast forward to mesh entry.\n");
    }
  }

  // Distance remaining to tally
  double distance_remaining = d_flight;

  while (distance_remaining > 0.) {
    // Distance we will travel in this cell
    // indices_tracklength.clear();
    auto next_tile = distance_to_next_index(r, u_, on, i, j, k);

    if (next_tile.first == INF) {
      // Something went wrong.... Don't score.
      Output::instance().save_warning(
          "Problem encountered with mesh tally in box_tally.\n");
      return indices_tracklength;
      // break;
    } else if (next_tile.first < 0.) {
      // Something went wrong.... Don't score.
      warning("Negative distance encountered with mesh tally.");
    }

    double d_tile = std::min(next_tile.first, distance_remaining);

    // Make the score if we are in a valid cell
    if (i >= 0 && i < static_cast<int>(Nx_) && j >= 0 &&
        j < static_cast<int>(Ny_) && k >= 0 && k < static_cast<int>(Nz_)) {
      size_t ui = static_cast<size_t>(i);
      size_t uj = static_cast<size_t>(j);
      size_t uk = static_cast<size_t>(k);

      StaticVector3 u_index = reduce_dimension(ui, uj, uk);
      trlen_d.index = u_index;
      trlen_d.distance = d_tile;
      indices_tracklength.push_back(trlen_d);

    } else {
      // If we arrive here, it means that we have left the tally region
      // when were we initially inside it. We can return here, as it's
      // impossible to go back in.
      return indices_tracklength;
    }

    // Remove the traveled distance
    distance_remaining -= d_tile;

    if (distance_remaining <= 0.) break;

    // Update the position and cell indices
    r = r + d_tile * u_;
    update_indices(next_tile.second, i, j, k, on);

  }  // While we still have to travel

  return indices_tracklength;
}

bool RegularCartesianMeshFilter::find_entry_point(Position& r,
                                                  const Direction& u,
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

void RegularCartesianMeshFilter::initialize_indices(
    const Position& r, const Direction& u, int& i, int& j, int& k,
    std::array<int, 3>& on) const {
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

void RegularCartesianMeshFilter::update_indices(int key, int& i, int& j, int& k,
                                                std::array<int, 3>& on) const {
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

std::pair<double, int> RegularCartesianMeshFilter::distance_to_next_index(
    const Position& r, const Direction& u, const std::array<int, 3>& on, int i,
    int j, int k) const {
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

void RegularCartesianMeshFilter::write_to_hdf5(H5::Group& grp) const {
  // Save id in attributes
  if (grp.hasAttribute("id")) {
    grp.deleteAttribute("id");
  }
  grp.createAttribute("id", this->id());

  // Save type in attributes
  if (grp.hasAttribute("type")) {
    grp.deleteAttribute("type");
  }
  grp.createAttribute("type", "regular-cartesian-mesh");

  // Save low position
  std::array<double, 3> r_low{r_low_.x(), r_low_.y(), r_low_.z()};
  if (grp.hasAttribute("low")) {
    grp.deleteAttribute("low");
  }
  grp.createAttribute("low", r_low);

  // Save high position
  std::array<double, 3> r_high{r_high_.x(), r_high_.y(), r_high_.z()};
  if (grp.hasAttribute("high")) {
    grp.deleteAttribute("high");
  }
  grp.createAttribute("high", r_high);

  // Save shape
  std::array<std::size_t, 3> shape{Nx_, Ny_, Nz_};
  if (grp.hasAttribute("shape")) {
    grp.deleteAttribute("shape");
  }
  grp.createAttribute("shape", shape);

  std::vector<double> x_bounds(Nx_ + 1, 0.);
  for (std::size_t i = 0; i <= Nx_; i++) {
    x_bounds[i] = (static_cast<double>(i) * dx_) + r_low_.x();
  }
  grp.createDataSet("x-bounds", x_bounds);

  std::vector<double> y_bounds(Ny_ + 1, 0.);
  for (std::size_t i = 0; i <= Ny_; i++) {
    y_bounds[i] = (static_cast<double>(i) * dy_) + r_low_.y();
  }
  grp.createDataSet("y-bounds", y_bounds);

  std::vector<double> z_bounds(Nz_ + 1, 0.);
  for (std::size_t i = 0; i <= Nz_; i++) {
    z_bounds[i] = (static_cast<double>(i) * dz_) + r_low_.z();
  }
  grp.createDataSet("z-bounds", z_bounds);
}

// Make the cartesian or position filter class
std::shared_ptr<RegularCartesianMeshFilter> make_regular_cartesian_mesh_filter(
    const YAML::Node& node) {
  if (!node["id"] || !node["id"].IsScalar()) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  if (!node["low"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"low\" coordinates are not provided.";
    fatal_error(mssg.str());
  } else if (!node["low"].IsSequence() || node["low"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", the given entry for the \"low\" coordinates must be a sequence "
            "of size 3.";
    fatal_error(mssg.str());
  }

  if (!node["high"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"high\" coordinates are not provided.";
    fatal_error(mssg.str());
  } else if (!node["high"].IsSequence() || node["high"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", the given entry for the \"high\" coordinates must be a sequence "
            "of size 3.";
    fatal_error(mssg.str());
  }

  if (!node["shape"]) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"shape\" is not provided.";
    fatal_error(mssg.str());
  } else if (!node["shape"].IsSequence() || node["shape"].size() != 3) {
    std::stringstream mssg;
    mssg << "For position-filter with id " << id
         << ", \"shape\" must be a sequence of size 3.";
    fatal_error(mssg.str());
  }

  std::vector<double> low_point = node["low"].as<std::vector<double>>();
  std::vector<double> high_point = node["high"].as<std::vector<double>>();
  std::vector<std::size_t> shape = node["shape"].as<std::vector<std::size_t>>();

  Position r_low(low_point[0], low_point[1], low_point[2]);
  Position r_high(high_point[0], high_point[1], high_point[2]);

  std::shared_ptr<RegularCartesianMeshFilter> mesh_type_filter =
      std::make_shared<RegularCartesianMeshFilter>(r_low, r_high, shape[0],
                                                   shape[1], shape[2], id);

  return mesh_type_filter;
}