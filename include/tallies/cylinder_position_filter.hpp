#ifndef CYLINDER_POSITION_FILTER_H
#define CYLINDER_POSITION_FILTER_H

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

#include <array>

class CylinderFilter : public PositionFilter {
 public:
  enum class Orientation { X, Y, Z };

  CylinderFilter(Position origin, double radius, double dx, double dy,
                 double dz, std::size_t nx, std::size_t ny, std::size_t nz,
                 Orientation z_, std::size_t id);

  StaticVector3 get_indices(const Tracker& tktr);

  StaticVector3 get_shape() {
    if ( Real_nx == 1 && Real_ny == 1 && Real_nz == 1){
      return {1};
    }
    StaticVector3 filter_shape{Nx_, Ny_, Nz_};

    map_indexes(filter_shape);
    return reduce_dimension(filter_shape[0], filter_shape[1], filter_shape[2]);
  }

  double radius() const { return radius_; }
  double inv_radius() const { return inv_radius_; }
  double get_scaled_radius(const Position& r) const;

  double z_min(const StaticVector3& indices) const;
  double z_max(const StaticVector3& indices) const;

  Position get_center(const StaticVector3& indices) const;

  std::pair<double, double> get_scaled_radius_and_angle(
      const StaticVector3& indices, const Position& r) const;

  std::string type_str() const override { return "cylinderpostionfilter"; }

  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& /*trkr*/, double /*d_flight*/) {
    std::vector<TracklengthDistance> index_;
    return index_;
  }

  Orientation get_axial_direction() { return length_axis_; }

 private:
  Position origin_, r_low_;
  std::size_t Nx_, Ny_, Nz_;
  std::size_t Real_nx, Real_ny, Real_nz;
  Orientation length_axis_;
  double radius_, pitch_x_, pitch_y_, dz_, inv_radius_, inv_pitch_x_,
      inv_pitch_y_, inv_dz_;

  // To map the indexes to either converting into class co-ordinate or into
  // original
  void map_indexes(StaticVector3& indexes) const {
    if (length_axis_ == Orientation::Z) {
      return;
    } else if (length_axis_ == Orientation::Y) {
      std::size_t nz = indexes[2];
      indexes[2] = indexes[1];
      indexes[1] = nz;
      return;
    } else if (length_axis_ == Orientation::X) {
      std::size_t nz = indexes[2];
      indexes[2] = indexes[0];
      indexes[0] = nz;
      return;
    }
  }
  // To map the positions to either converting into class co-ordinate or into
  // original
  Position map_coordinate(const Position& point) const {
    if (length_axis_ == Orientation::Z) {
      return point;
    } else if (length_axis_ == Orientation::Y) {
      return Position(point.x(), point.z(), point.y());
    } else if (length_axis_ == Orientation::X) {
      return Position(point.z(), point.y(), point.x());
    }
    return point;
  }

  // function will reduce the dimsion, if there is only one bin in the direction
  StaticVector3 reduce_dimension(const size_t& loc_x, const size_t& loc_y,
                                 const size_t& loc_z) {
    StaticVector3 reduce_;
    if ( Real_nx == 1 && Real_ny == 1 && Real_nz == 1){
      return {loc_x};
    }
    if (Real_nx > 1) {
      reduce_.push_back(loc_x);
    }

    if (Real_ny > 1) {
      reduce_.push_back(loc_y);
    }

    if (Real_nz > 1) {
      reduce_.push_back(loc_z);
    }
    return reduce_;
  }
};

std::shared_ptr<CylinderFilter> make_cylinder_position_filter(
    const YAML::Node& node);

#endif