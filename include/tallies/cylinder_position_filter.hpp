#ifndef CYLINDER_POSITION_FILTER_H
#define CYLINDER_POSITION_FILTER_H

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

#include <array>

class CylinderPositionFilter : public PositionFilter {
 public:
  enum class Orientation { X, Y, Z };

  CylinderPositionFilter(Position origin, double radius, double dx, double dy,
                         double dz, std::size_t nx, std::size_t ny,
                         std::size_t nz, Orientation z_, std::size_t id);

  StaticVector3 get_indices(const Tracker& tktr);

  StaticVector3 get_shape() {
    StaticVector3 filter_shape{Nx_, Ny_, Nz_};

    map_indexes(filter_shape);
    return filter_shape;
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

 private:
  Position origin_, r_low_;
  std::size_t Nx_, Ny_, Nz_;
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
};

std::shared_ptr<CylinderPositionFilter> make_cylinder_position_filter(
    const YAML::Node& node);

#endif