#include <tallies/cylinder_position_filter.hpp>
#include <utils/error.hpp>

#include <vector>

CylinderPositionFilter::CylinderPositionFilter(Position origin, double radius,
                                               double dx, double dy, double dz,
                                               std::size_t nx, std::size_t ny,
                                               std::size_t nz, Orientation z_,
                                               std::size_t id)
    : PositionFilter(id),
      origin_(origin),
      r_low_(),
      Nx_(nx),
      Ny_(ny),
      Nz_(nz),
      length_axis_(z_),
      radius_(radius),
      pitch_x_(dx),
      pitch_y_(dy),
      dz_(dz),
      inv_radius_(),
      inv_pitch_x_(),
      inv_pitch_y_(),
      inv_dz_() {
  // Map the parameters according to the orientation
  // since this class can do the caluclation assumind z-axis as axial
  // therefore, all the parameter will be mapped to z-axis,
  // however, while return out from this class, parameters
  // will be mapped back to its original orientation
  origin_ = map_coordinate(origin);
  std::string str_axial = "z";
  std::string str_pitch_x = "x";
  std::string str_pitch_y = "y";
  if (length_axis_ == Orientation::X) {
    pitch_x_ = dz;
    dz_ = dx;
    Nx_ = nz;
    Nz_ = nx;
    str_axial = "x";
    str_pitch_x = "z";
  }

  if (length_axis_ == Orientation::Y) {
    pitch_y_ = dz;
    dz_ = dy;
    Nz_ = ny;
    Ny_ = nz;
    str_axial = "y";
    str_pitch_y = "z";
  }

  // Check for the valid radius
  if (radius_ <= 0.)
    fatal_error("A non-zero value of radius must be provided.");
  inv_radius_ = 1.0 / radius_;

  // Check for a non-zero pitch_x_ and pitch_y_ should be more than cylinder
  // diameter
  if ((pitch_x_ == 0.) && (Nx_ != 1)) {
    warning(str_pitch_x +
            "-direction pitch is zero, Only one will be considered in "
            "that-direction.");

    Nx_ = static_cast<size_t>(1);
    inv_pitch_x_ = 0.;
  } else if (pitch_x_ != 0. && Nx_ == 1) {
    warning(str_pitch_x +
            "-direction pitch is non-zero, but the bin is one in that "
            "direction,so, Only one will be considered in "
            "that-direction.");
  } else {
    if (pitch_x_ != 0. && pitch_x_ < (2. * radius_))
      fatal_error(str_pitch_x + "-pitch is less than the radius");
    inv_pitch_x_ = 1.0 / pitch_x_;
  }

  if ((pitch_y_ == 0.) && (Ny_ != 1)) {
    warning(str_pitch_y +
            "-direction pitch is zero, Only one will be considered in "
            "the direction.");

    Ny_ = static_cast<size_t>(1);
    inv_pitch_y_ = 0.;
  } else if (pitch_y_ != 0. && Ny_ == 1) {
    warning(str_pitch_y +
            "-direction pitch is non-zero, but the bin is one in that "
            "direction,so, Only one will be considered in "
            "that-direction.");
  } else {
    if (pitch_y_ != 0. && pitch_y_ < (2. * radius_))
      fatal_error(str_pitch_y + "pitch is less than the radius ");
    inv_pitch_y_ = 1.0 / pitch_y_;
  }

  // A finite length the z-direction should be given.
  if (dz_ == 0.)
    fatal_error("Length of the cylinder shouldn't zero.");
  else
    inv_dz_ = 1.0 / dz_;

  // calculate low point for the purpose of getting the indices
  // r_low_ is mapped, means z-cooridnate will always be axial direction
  double low_x = origin_.x() - pitch_x_ * 0.5;
  double low_y = origin_.y() - pitch_y_ * 0.5;
  r_low_ = Position(low_x, low_y, origin_.z());
}

StaticVector3 CylinderPositionFilter::get_indices(const Tracker& tktr) {
  const Position r = map_coordinate(tktr.r());
  int nx = 0, ny = 0;

  if (pitch_x_ != 0.)
    nx = static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_pitch_x_));

  if (pitch_y_ != 0.)
    ny = static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_pitch_y_));

  int nz = static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_dz_));

  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(nx);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(ny);

  StaticVector3 indices;
  // check if nx, ny, and nz are positive
  // and less the number of bins in that direction
  if (nx >= 0 && nx < static_cast<int>(Nx_) && ny >= 0 &&
      ny < static_cast<int>(Ny_) && nz >= 0 && nz < static_cast<int>(Nz_)) {
    // check if the position is inside the circular radius or not
    if (sqrt((new_origin_x - r.x()) * (new_origin_x - r.x()) +
             (new_origin_y - r.y()) * (new_origin_y - r.y())) <=
        (radius_ + 1E-15)) {
      indices.push_back(static_cast<std::size_t>(nx));
      indices.push_back(static_cast<std::size_t>(ny));
      indices.push_back(static_cast<std::size_t>(nz));
      map_indexes(indices);
    } else {
      return indices;
    }
  } else {
    return indices;
  }

  return indices;
}

double CylinderPositionFilter::z_min(const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  return (r_low_.z() + static_cast<double>(index[2]) * dz_);
}

double CylinderPositionFilter::z_max(const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  return (r_low_.z() + static_cast<double>(index[2]) * dz_ + dz_);
}

Position CylinderPositionFilter::get_center(
    const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(index[0]);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(index[1]);
  double new_origin_z = origin_.z() + dz_ * static_cast<double>(index[2]);
  return map_coordinate(Position(new_origin_x, new_origin_y, new_origin_z));
}

// first will be the scaled radius and second will be the angle
std::pair<double, double> CylinderPositionFilter::get_scaled_radius_and_angle(
    const StaticVector3& indices, const Position& r) const {
  // the indicies and Position is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  Position mapped_r = map_coordinate(r);
  const double rx = r.x();
  const double ry = r.y();
  // get new centre points
  const double new_origin_x =
      origin_.x() + pitch_x_ * static_cast<double>(index[0]);
  const double new_origin_y =
      origin_.y() + pitch_y_ * static_cast<double>(index[1]);
  // scaled-radius
  const double scaled_r = (std::sqrt(rx * rx + ry * ry)) * inv_radius_;
  // for tan(theta) = x/ y; get x and y
  const double base = rx - new_origin_x;
  const double height = ry - new_origin_y;
  const double theta = atan(height / base);
  return {scaled_r, theta};
}

std::shared_ptr<CylinderPositionFilter> make_cylinder_position_filter(
    const YAML::Node& node) {
  // get the id
  if (!node["id"] || !node["id"].IsScalar()) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  // check and get the origin
  if (!node["origin"]) {
    fatal_error(
        "the origin must be provided for the cylinder position filter.");
  } else if (!node["origin"].IsSequence() || node["origin"].size() != 3) {
    fatal_error("Invalid origin coordinated are given.");
  }
  std::vector<double> origin_point = node["origin"].as<std::vector<double>>();
  Position origin(origin_point[0], origin_point[1], origin_point[2]);

  // Get the radius
  if (!node["radius"]) {
    fatal_error("the radius must be given for the cylinder position filter.");
  } else if (!node["radius"].IsScalar()) {
    fatal_error("Given radius must be a scalar only.");
  }
  double radius = node["radius"].as<double>();

  // get the length of the cylinder
  if (!node["axial-length"]) {
    fatal_error("the axial length must be given.");
  } else if (!node["axial-length"].IsScalar()) {
    fatal_error("Invalid axial-length is given.");
  }
  double length = node["axial-length"].as<double>();

  // get the axial-direction axis name
  if (!node["axial-axis"]) {
    fatal_error(
        "axial-direction axis must be given for the cylinder position filter.");
  } else if (!node["axial-axis"].IsScalar()) {
    fatal_error(
        "Invalid axial-axis is given for the cylinder position filter.");
  }
  std::string axial_axis = node["axial-axis"].as<std::string>();

  // the default shape will nx = 1, ny = 1, and nz = 1;
  std::size_t nx = 1;
  std::size_t ny = 1;
  std::size_t nz = 1;
  // get the shape
  if (!node["shape"]) {
    warning("the shape must be given for the cylinder position filter");
  } else if (!node["shape"].IsSequence() || node["shape"].size() != 3) {
    fatal_error("Invalid shape is given for the cylinder position filter");
  } else {
    std::vector<std::size_t> shape =
        node["shape"].as<std::vector<std::size_t>>();
    nx = shape[0];
    ny = shape[1];
    nz = shape[2];
  }

  // get the pitches for each
  std::string str_pitch_x, str_pitch_y;
  if (axial_axis == "z" || axial_axis == "Z") {
    str_pitch_x = "pitch-x";
    str_pitch_y = "pitch-y";
  } else if (axial_axis == "x" || axial_axis == "X") {
    str_pitch_x = "pitch-z";
    str_pitch_y = "pitch-y";
  } else if (axial_axis == "y" || axial_axis == "Y") {
    str_pitch_x = "pitch-x";
    str_pitch_y = "pitch-z";
  }

  // get the pitch-x
  double pitch_x = 0.;
  if (!node[str_pitch_x]) {
    if (nx != 1) {
      warning(str_pitch_x +
              " is not given for the cylinder position filter, so there will "
              "be one cylinder is x-direction.");
    }
    nx = 1;
  } else {
    if (!node[str_pitch_x].IsScalar()) {
      fatal_error("Invalid " + str_pitch_x +
                  " is given for the cylinder position filter.");
    }
    pitch_x = node[str_pitch_x].as<double>();
    if (pitch_x == 0. && nx != 1) {
      warning("for the cylinder position filter, the " + str_pitch_x +
              " is 0.0, so, there will be 1 cylinder in x-direction.");
      nx = 1;
    }
  }

  // get the pitch-y
  double pitch_y = 0.;
  if (!node[str_pitch_y]) {
    if (nx != 1) {
      warning(str_pitch_y +
              " is not given for the cylinder position filter, so there will "
              "be one cylinder is x-direction.");
    }
    nx = 1;
  } else {
    if (!node[str_pitch_y].IsScalar()) {
      fatal_error("Invalid " + str_pitch_y +
                  " is given for the cylinder position filter.");
    }
    pitch_y = node[str_pitch_y].as<double>();
    if (pitch_y == 0. && nx != 1) {
      warning("for the cylinder position filter, the " + str_pitch_y +
              " is 0.0, so, there will be 1 cylinder in x-direction.");
      nx = 1;
    }
  }

  // make the cylinder position filter
  double dx, dy, dz;
  std::shared_ptr<CylinderPositionFilter> cylinder_filter;
  if (axial_axis == "x" || axial_axis == "X") {
    dx = length;
    dy = pitch_y;
    dz = pitch_x;
    cylinder_filter = std::make_shared<CylinderPositionFilter>(
        origin, radius, dx, dy, dz, nx, ny, nz,
        CylinderPositionFilter::Orientation::X, id);
  } else if (axial_axis == "y" || axial_axis == "Y") {
    dx = pitch_x;
    dy = length;
    dz = pitch_y;
    cylinder_filter = std::make_shared<CylinderPositionFilter>(
        origin, radius, dx, dy, dz, nx, ny, nz,
        CylinderPositionFilter::Orientation::Y, id);
  } else if (axial_axis == "z" || axial_axis == "Z") {
    dx = pitch_x;
    dy = pitch_y;
    dz = length;
    cylinder_filter = std::make_shared<CylinderPositionFilter>(
        origin, radius, dx, dy, dz, nx, ny, nz,
        CylinderPositionFilter::Orientation::Z, id);
  } else {
    fatal_error("Invalid axial direction axis is given.");
  }

  return cylinder_filter;
}