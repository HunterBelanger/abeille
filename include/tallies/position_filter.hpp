#ifndef POSITION_FILTER_H
#define POSITION_FILTER_H

#include <memory>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>
#include <utils/position.hpp>

#include <boost/container/static_vector.hpp>


using StaticVector6 = boost::container::static_vector<size_t, 6>;


enum class FilterType {
  Energy_Filter,
  Position_Filter,
  Cartesian_Filter,
  Box_Position_Filter,
  Mesh_Positin_Filter,
  Cylinder_Position_Filter,
  Cylinder_Array_Filter
};

struct TracklengthDistance {
  std::vector<size_t> indexes_;
  double distance_in_bin;
};

class PositionFilter {
 public:
  PositionFilter() = default;

  virtual ~PositionFilter() = default;

  // virtual bool get_indices(const Tracker& tktr, std::array<int, 3>& indices)
  // = 0;

  //virtual std::vector<size_t> get_indices(const Tracker& tktr) = 0;

  virtual StaticVector6 get_indices(const Tracker& tktr) = 0;

  virtual std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) = 0;

  virtual size_t Nx() const = 0;
  virtual size_t Ny() const = 0;
  virtual size_t Nz() const = 0;

  virtual StaticVector6 get_dimension() = 0;

  // Perhaps Not Needed
  virtual FilterType type() const { return FilterType::Position_Filter; }
  virtual std::string type_str() const = 0;
};

std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node);

#endif