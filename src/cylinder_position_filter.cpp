#include <tallies/cylinder_position_filter.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <vector>

CylinderFilter::CylinderFilter(Position origin, double radius, double dx,
                               double dy, double dz, std::size_t nx,
                               std::size_t ny, std::size_t nz, Orientation z_,
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
  } else if (length_axis_ == Orientation::Y) {
    pitch_y_ = dz;
    dz_ = dy;
    Nz_ = ny;
    Ny_ = nz;
    str_axial = "y";
    str_pitch_y = "z";
  }
  std::string id_str = std::to_string(id);
  // Check for the valid radius
  if (radius_ <= 0.)
    fatal_error("A non-zero value of radius must be provided for the id: " +
                id_str);
  inv_radius_ = 1.0 / radius_;

  // Check for a non-zero pitch_x_ and pitch_y_ should be more than cylinder
  // diameter
  if ((pitch_x_ == 0.) && (Nx_ != 1)) {
    std::stringstream messg;
    messg << "for id: " << id << ", " << str_pitch_x
          << "-direction pitch is zero, therefore, only one cylinder will be "
             "considered in "
          << str_pitch_x << "-direction";
    warning(messg.str());

    Nx_ = static_cast<size_t>(1);
    inv_pitch_x_ = 0.;
  } else if (pitch_x_ != 0. && Nx_ == 1) {
    std::stringstream messg;
    messg << "for id: " << id << ", " << str_pitch_x
          << "-direction pitch is a non-zero value, but the bin size is one, "
             "so, one cylinder will be considered in "
          << str_pitch_x << "-direction";
    warning(messg.str());
  } else {
    if (pitch_x_ != 0. && pitch_x_ < (2. * radius_))
      fatal_error("for id: " + id_str + ", " + str_pitch_x +
                  "-pitch is less than the radius");
    inv_pitch_x_ = 1.0 / pitch_x_;
  }

  if ((pitch_y_ == 0.) && (Ny_ != 1)) {
    std::stringstream messg;
    messg << "for id: " << id << ", " << str_pitch_y
          << "-direction pitch is zero, therefore, only one cylinder will be "
             "considered in "
          << str_pitch_y << "-direction";
    warning(messg.str());

    Ny_ = static_cast<size_t>(1);
    inv_pitch_y_ = 0.;
  } else if (pitch_y_ != 0. && Ny_ == 1) {
    std::stringstream messg;
    messg << "for id: " << id << ", " << str_pitch_y
          << "-direction pitch is a non-zero value, but the bin size is one, "
             "so, one cylinder will be considered in "
          << str_pitch_y << "-direction";
    warning(messg.str());
  } else {
    if (pitch_y_ != 0. && pitch_y_ < (2. * radius_))
      fatal_error("for id: " + id_str + ", " + str_pitch_y +
                  "pitch is less than the radius ");
    inv_pitch_y_ = 1.0 / pitch_y_;
  }

  // A finite length the z-direction should be given.
  if (dz_ == 0.)
    fatal_error("for id: " + id_str +
                ",length of the cylinder shouldn't zero.");
  else
    inv_dz_ = 1.0 / dz_;

  // calculate low point for the purpose of getting the indices
  // r_low_ is mapped, means z-cooridnate will always be axial direction
  double low_x = origin_.x() - pitch_x_ * 0.5;
  double low_y = origin_.y() - pitch_y_ * 0.5;
  r_low_ = Position(low_x, low_y, origin_.z());
}

StaticVector3 CylinderFilter::get_indices(const Tracker& tktr) {
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
      return indices;
    }
  }
  return indices;
}

double CylinderFilter::z_min(const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  return (r_low_.z() + static_cast<double>(index[2]) * dz_);
}

double CylinderFilter::z_max(const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  return (r_low_.z() + static_cast<double>(index[2]) * dz_ + dz_);
}

Position CylinderFilter::get_center(const StaticVector3& indices) const {
  // the indicies is not according to the class, so map it.
  StaticVector3 index = indices;
  map_indexes(index);
  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(index[0]);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(index[1]);
  double new_origin_z = origin_.z() + dz_ * static_cast<double>(index[2]);
  return map_coordinate(Position(new_origin_x, new_origin_y, new_origin_z));
}

// first will be the scaled radius and second will be the angle
std::pair<double, double> CylinderFilter::get_scaled_radius_and_angle(
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
  const double theta = std::atan(height / base);
  return {scaled_r, theta};
}

std::shared_ptr<CylinderFilter> make_cylinder_position_filter(
    const YAML::Node& node) {
  // get the id
  if (!node["id"] || !node["id"].IsScalar()) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();
  std::string id_str = std::to_string(id);

  // check and get the origin
  if (!node["origin"]) {
    fatal_error(
        "for id: " + id_str +
        "the origin must be provided for the cylinder position filter.");
  } else if (!node["origin"].IsSequence() || node["origin"].size() != 3) {
    fatal_error("for id: " + id_str +
                ", invalid origin coordinated are given.");
  }
  std::vector<double> origin_point = node["origin"].as<std::vector<double>>();
  Position origin(origin_point[0], origin_point[1], origin_point[2]);

  // Get the radius
  if (!node["radius"]) {
    fatal_error("for id: " + id_str +
                ", the radius must be given for the cylinder position filter.");
  } else if (!node["radius"].IsScalar()) {
    fatal_error("for id: " + id_str + ", given radius must be a scalar only.");
  }
  double radius = node["radius"].as<double>();

  // get the length of the cylinder
  if (!node["length"]) {
    fatal_error("for id: " + id_str + ", the axial length must be given.");
  } else if (!node["length"].IsScalar()) {
    fatal_error("for id: " + id_str + ", invalid length is given.");
  }
  double length = node["length"].as<double>();

  // get the axial-direction axis name
  if (!node["axis"]) {
    fatal_error(
        "for id: " + id_str +
        "axial-direction axis must be given for the cylinder position filter.");
  } else if (!node["axis"].IsScalar()) {
    fatal_error("for id: " + id_str +
                ", invalid axis is given for the cylinder position filter.");
  }
  std::string axial_axis = node["axis"].as<std::string>();

  // the default shape will nx = 1, ny = 1, and nz = 1;
  std::size_t nx = 1;
  std::size_t ny = 1;
  std::size_t nz = 1;
  // get the shape
  if (!node["shape"]) {
    warning("for id: " + id_str +
            ", the shape must be given for the cylinder position filter");
  } else if (!node["shape"].IsSequence() || node["shape"].size() != 3) {
    fatal_error("for id: " + id_str +
                ", invalid shape is given for the cylinder position filter");
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
              "be one cylinder is x-direction, for id: " +
              id_str + "");
    }
    nx = 1;
  } else {
    if (!node[str_pitch_x].IsScalar()) {
      fatal_error("for id: " + id_str + ", invalid " + str_pitch_x +
                  " is given for the cylinder position filter.");
    }
    pitch_x = node[str_pitch_x].as<double>();
    if (pitch_x == 0. && nx != 1) {
      warning("for id: " + id_str + ", given" + str_pitch_x +
              " is 0.0, so, there will be 1 cylinder in x-direction.");
      nx = 1;
    }
  }

  // get the pitch-y
  double pitch_y = 0.;
  if (!node[str_pitch_y]) {
    if (ny != 1) {
      warning(str_pitch_y +
              " is not given for the cylinder position filter, so there will "
              "be one cylinder is x-direction, for id: " +
              id_str + "");
    }
    ny = 1;
  } else {
    if (!node[str_pitch_y].IsScalar()) {
      fatal_error("for id: " + id_str + ", invalid " + str_pitch_y +
                  " is given for the cylinder position filter.");
    }
    pitch_y = node[str_pitch_y].as<double>();
    if (pitch_y == 0. && ny != 1) {
      warning("for id: " + id_str + ", given" + str_pitch_x +
              " is 0.0, so, there will be 1 cylinder in x-direction.");
      ny = 1;
    }
  }

  // make the cylinder position filter
  double dx, dy, dz;
  std::shared_ptr<CylinderFilter> cylinder_filter;
  if (axial_axis == "x" || axial_axis == "X") {
    dx = length;
    dy = pitch_y;
    dz = pitch_x;
    cylinder_filter =
        std::make_shared<CylinderFilter>(origin, radius, dx, dy, dz, nx, ny, nz,
                                         CylinderFilter::Orientation::X, id);
  } else if (axial_axis == "y" || axial_axis == "Y") {
    dx = pitch_x;
    dy = length;
    dz = pitch_y;
    cylinder_filter =
        std::make_shared<CylinderFilter>(origin, radius, dx, dy, dz, nx, ny, nz,
                                         CylinderFilter::Orientation::Y, id);
  } else if (axial_axis == "z" || axial_axis == "Z") {
    dx = pitch_x;
    dy = pitch_y;
    dz = length;
    cylinder_filter =
        std::make_shared<CylinderFilter>(origin, radius, dx, dy, dz, nx, ny, nz,
                                         CylinderFilter::Orientation::Z, id);
  } else {
    fatal_error("For id: " + id_str +
                ", invalid axial direction axis is given.");
  }

  return cylinder_filter;
}