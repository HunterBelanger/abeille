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

  StaticVector3 get_indices(const Tracker& tktr) const override final;

  StaticVector3 get_position_index(const Position& r) const;

  StaticVector3 get_shape() const override final;

  double radius() const { return radius_; }
  double inv_radius() const { return inv_radius_; }
  double get_scaled_radius(const Position& r) const;

  double z_min(const StaticVector3& indices) const;
  double z_max(const StaticVector3& indices) const;
  double dz() const { return dz_; }
  double inv_dz() const { return inv_dz_; }
  bool is_infinite_cylinder() const { return infinite_length_; }

  Position get_center(const StaticVector3& indices, const bool is_map) const;

  std::pair<double, double> get_scaled_radius_and_angle(
      const StaticVector3& indices, const Position& r) const;

  std::string type_str() const override { return "cylinder-filter"; }

  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) const override final;

  Orientation get_axial_direction() { return length_axis_; }

  void write_to_hdf5(H5::Group& grp) const override final;

 private:
  Position origin_, r_low_;
  std::size_t Nx_, Ny_, Nz_;
  std::size_t Real_nx_, Real_ny_, Real_nz_;
  Orientation length_axis_;
  double radius_, pitch_x_, pitch_y_, dz_, inv_radius_, inv_pitch_x_,
      inv_pitch_y_, inv_dz_;
  std::size_t x_index, y_index,
      z_index;  // to access the correct location of the indices
  // these index-locations are based on the real-position of x, y, z. Not based
  // on the class orientation
  bool infinite_length_ = false;  // to incorporate the infinte cylinder
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
                                 const size_t& loc_z) const {
    StaticVector3 reduce_;
    if (Real_nx_ == 1 && Real_ny_ == 1 && Real_nz_ == 1) {
      return {1};
    }

    if (infinite_length_ == true) {
      if (Nx_ == 1 && Ny_ == 1) {
        return {1};
      }
    }

    if (Real_nx_ > 1) {
      reduce_.push_back(loc_x);
    }

    if (Real_ny_ > 1) {
      reduce_.push_back(loc_y);
    }

    if (Real_nz_ > 1) {
      reduce_.push_back(loc_z);
    }
    return reduce_;
  }

};

std::shared_ptr<CylinderFilter> make_cylinder_filter(const YAML::Node& node);

#endif