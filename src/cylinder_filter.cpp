#include <tallies/cylinder_filter.hpp>
#include <utils/error.hpp>

#include <cmath>
#include <sstream>
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
      Real_nx(nx),
      Real_ny(ny),
      Real_nz(nz),
      length_axis_(z_),
      radius_(radius),
      pitch_x_(dx),
      pitch_y_(dy),
      dz_(dz),
      inv_radius_(1. / radius_),
      inv_pitch_x_(),
      inv_pitch_y_(),
      inv_dz_() {
  // Map the parameters according to the orientation
  // since this class can do the caluclation assumind z-axis as axial
  // therefore, all the parameter will be mapped to z-axis,
  // however, while return out from this class, parameters
  // will be mapped back to its original orientation
  origin_ = map_coordinate(origin);
  if (length_axis_ == Orientation::X) {
    pitch_x_ = dz;
    dz_ = dx;
    Nx_ = nz;
    Nz_ = nx;
  } else if (length_axis_ == Orientation::Y) {
    pitch_y_ = dz;
    dz_ = dy;
    Nz_ = ny;
    Ny_ = nz;
  }

  // Check for the valid radius
  if (radius_ <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " was provided with a negative or zero radius.";
    fatal_error(mssg.str());
  }

  // Make sure pitches are all >= 0.
  if (dx <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dx <= 0.";
    fatal_error(mssg.str());
  }

  if (dy <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dy <= 0.";
    fatal_error(mssg.str());
  }

  if (dz <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with a dz <= 0.";
    fatal_error(mssg.str());
  }

  // Make sure shapes are all > 0
  if (nx == 0 && length_axis_ == Orientation::X) {
    infinite_length = true;
  } else if (nx == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with nx = 0.";
    fatal_error(mssg.str());
  }

  if (ny == 0 && length_axis_ == Orientation::Y) {
    infinite_length = true;
  } else if (ny == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with ny = 0.";
    fatal_error(mssg.str());
  }

  if (nz == 0 && length_axis_ == Orientation::Z) {
    infinite_length = true;
  } else if (nz == 0) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " provided with nz = 0.";
    fatal_error(mssg.str());
  }

  // Make sure the pitch_x and pitch_y are >= diameter
  if (pitch_x_ < 2. * radius_ || pitch_y_ < 2. * radius_) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " has dimensions which are shorter than the diameter.";
    fatal_error(mssg.str());
  }

  // Calculate inverse values
  inv_pitch_x_ = 1.0 / pitch_x_;
  inv_pitch_y_ = 1.0 / pitch_y_;
  inv_dz_ = 1.0 / dz_;

  // calculate low point for the purpose of getting the indices
  // r_low_ is mapped, means z-cooridnate will always be axial direction
  double low_x = origin_.x() - pitch_x_ * 0.5;
  double low_y = origin_.y() - pitch_y_ * 0.5;
  r_low_ = Position(low_x, low_y, origin_.z());
}

StaticVector3 CylinderFilter::get_indices(const Tracker& tktr) const {
  const Position r = map_coordinate(tktr.r());

  const int nx =
      static_cast<int>(std::floor((r.x() - r_low_.x()) * inv_pitch_x_));
  const int ny =
      static_cast<int>(std::floor((r.y() - r_low_.y()) * inv_pitch_y_));
  const int nz =
      infinite_length
          ? 0
          : static_cast<int>(std::floor((r.z() - r_low_.z()) * inv_dz_));

  double new_origin_x = origin_.x() + pitch_x_ * static_cast<double>(nx);
  double new_origin_y = origin_.y() + pitch_y_ * static_cast<double>(ny);

  StaticVector3 indices;
  // check if nx, ny, and nz are positive
  // and less the number of bins in that direction
  if (nx >= 0 && nx < static_cast<int>(Nx_) && ny >= 0 &&
      ny < static_cast<int>(Ny_) &&
      ((nz >= 0 && nz < static_cast<int>(Nz_)) || infinite_length)) {
    // check if the position is inside the circular radius or not
    if (sqrt((new_origin_x - r.x()) * (new_origin_x - r.x()) +
             (new_origin_y - r.y()) * (new_origin_y - r.y())) <=
        (radius_ + 1E-15)) {
      indices.push_back(static_cast<std::size_t>(nx));
      indices.push_back(static_cast<std::size_t>(ny));
      indices.push_back(static_cast<std::size_t>(nz));
      map_indexes(indices);
      return reduce_dimension(indices[0], indices[1], indices[2]);
    }
  }
  return indices;
}

StaticVector3 CylinderFilter::get_shape() const {
  if (Real_nx == 1 && Real_ny == 1 && Real_nz == 1) {
    return {1};
  }
  StaticVector3 filter_shape{Nx_, Ny_, Nz_};

  map_indexes(filter_shape);
  return reduce_dimension(filter_shape[0], filter_shape[1], filter_shape[2]);
}

std::vector<TracklengthDistance> CylinderFilter::get_indices_tracklength(
    const Tracker& /*trkr*/, double /*d_flight*/) const {
  fatal_error("Not yet implemented.");

  return {};
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

void CylinderFilter::write_to_hdf5(H5::Group& grp) const {
  // Save id in attributes
  if (grp.hasAttribute("id")) {
    grp.deleteAttribute("id");
  }
  grp.createAttribute("id", this->id());

  // Save type in attributes
  if (grp.hasAttribute("type")) {
    grp.deleteAttribute("type");
  }
  grp.createAttribute("type", "cylinder-filter");

  // Save origin position
  std::array<double, 3> origin{origin_.x(), origin_.y(), origin_.z()};
  if (grp.hasAttribute("origin")) {
    grp.deleteAttribute("origin");
  }
  grp.createAttribute("origin", origin);

  // Save radius
  if (grp.hasAttribute("radius")) {
    grp.deleteAttribute("radius");
  }
  grp.createAttribute("radius", this->radius());

  // Save axis
  if (grp.hasAttribute("axis")) {
    grp.deleteAttribute("axis");
  }
  switch (length_axis_) {
    case Orientation::X:
      grp.createAttribute("axis", "x");
      break;

    case Orientation::Y:
      grp.createAttribute("axis", "y");
      break;

    case Orientation::Z:
      grp.createAttribute("axis", "z");
      break;
  }

  // Save shape
  std::array<std::size_t, 3> shape{Real_nx, Real_ny, Real_nz};
  switch (length_axis_) {
    case Orientation::X:
      grp.createAttribute("axis", "x");
      break;

    case Orientation::Y:
      grp.createAttribute("axis", "y");
      break;

    case Orientation::Z:
      grp.createAttribute("axis", "z");
      break;
  }
  if (grp.hasAttribute("shape")) {
    grp.deleteAttribute("shape");
  }
  grp.createAttribute("shape", shape);

  // Save pitch
  std::array<double, 3> pitch{pitch_x_, pitch_y_, dz_};
  if (length_axis_ == Orientation::X) {
    pitch[0] = dz_;
    pitch[2] = pitch_x_;
  } else if (length_axis_ == Orientation::Y) {
    pitch[1] = dz_;
    pitch[2] = pitch_y_;
  }
  if (grp.hasAttribute("pitch")) {
    grp.deleteAttribute("pitch");
  }
  grp.createAttribute("pitch", pitch);
}

std::shared_ptr<CylinderFilter> make_cylinder_filter(const YAML::Node& node) {
  // Get the id
  if (!node["id"] || node["id"].IsScalar() == false) {
    fatal_error("Invalid id is given for the position-filter.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  // check and get the origin
  if (!node["origin"] || node["origin"].IsSequence() == false ||
      node["origin"].size() != 3) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid origin entry.";
    fatal_error(mssg.str());
  }
  std::vector<double> origin_point = node["origin"].as<std::vector<double>>();
  Position origin(origin_point[0], origin_point[1], origin_point[2]);

  // Get the radius
  if (!node["radius"] || node["radius"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid radius entry.";
    fatal_error(mssg.str());
  }
  const double radius = node["radius"].as<double>();
  if (radius <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has a radius which is <= 0.";
    fatal_error(mssg.str());
  }

  // Get the axial-direction axis name
  if (!node["axis"] || node["axis"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " is missing a valid axis entry.";
    fatal_error(mssg.str());
  }
  const std::string axial_axis = node["axis"].as<std::string>();
  if (axial_axis != "x" && axial_axis != "y" && axial_axis != "z") {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " was provided with invalid axis entry \"" << axial_axis << "\".";
    fatal_error(mssg.str());
  }
  CylinderFilter::Orientation orientation;
  if (axial_axis == "x")
    orientation = CylinderFilter::Orientation::X;
  else if (axial_axis == "y")
    orientation = CylinderFilter::Orientation::Y;
  else
    orientation = CylinderFilter::Orientation::Z;

  // The default shape will nx = 1, ny = 1, and nz = 1;
  std::size_t nx = 1;
  std::size_t ny = 1;
  std::size_t nz = 1;
  // Get the shape
  if (node["shape"] && node["shape"].IsSequence() &&
      node["shape"].size() == 3) {
    std::vector<std::size_t> shape =
        node["shape"].as<std::vector<std::size_t>>();
    nx = shape[0];
    ny = shape[1];
    nz = shape[2];
  } else if (node["shape"]) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has an invalid shape entry.";
    fatal_error(mssg.str());
  }

  // Get the pitch of the lattice/box.
  if (!node["pitch"] || node["pitch"].IsSequence() == false ||
      node["pitch"].size() != 3) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has invalid pitch entry.";
    fatal_error(mssg.str());
  }
  std::vector<double> pitches = node["pitch"].as<std::vector<double>>();

  double pitch_x, pitch_y, length;
  if (axial_axis == "z") {
    pitch_x = pitches[0];
    pitch_y = pitches[1];
    length = pitches[2];
  } else if (axial_axis == "x") {
    pitch_x = pitches[2];
    pitch_y = pitches[1];
    length = pitches[0];
  } else {
    pitch_x = pitches[0];
    pitch_y = pitches[2];
    length = pitches[1];
  }

  if (pitch_x <= 0. || pitch_y <= 0. || length <= 0.) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id << " has pitches which are <= 0.";
    fatal_error(mssg.str());
  }

  if (pitch_x < 2. * radius || pitch_y < 2. * radius) {
    std::stringstream mssg;
    mssg << "Cylinder filter with id " << id
         << " has pitches which are too small for the given radius.";
    fatal_error(mssg.str());
  }

  // Make the filter
  return std::make_shared<CylinderFilter>(origin, radius, pitches[0],
                                          pitches[1], pitches[2], nx, ny, nz,
                                          orientation, id);
}