#ifndef POSITION_FILTER_H
#define POSITION_FILTER_H

#include <simulation/tracker.hpp>



#include <yaml-cpp/yaml.h>
#include <boost/container/static_vector.hpp>

#include <memory>
#include <vector>

using StaticVector3 = boost::container::static_vector<size_t, 3>;

struct TracklengthDistance {
  std::vector<size_t> indexes_;
  double distance_in_bin;
};

class PositionFilter {
 public:
  PositionFilter() = default;

  virtual ~PositionFilter() = default;

  virtual StaticVector3 get_indices(const Tracker& tktr) = 0;

  virtual std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) = 0;

  virtual size_t Nx() const = 0;
  virtual size_t Ny() const = 0;
  virtual size_t Nz() const = 0;

  virtual StaticVector3 get_dimension() = 0;

  // Perhaps Not Needed

  virtual std::string type_str() const = 0;
};

std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node);

#endif