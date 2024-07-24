#ifndef POSITION_FILTER_H
#define POSITION_FILTER_H

#include <simulation/tracker.hpp>

#include <yaml-cpp/yaml.h>
#include <boost/container/static_vector.hpp>
using StaticVector3 = boost::container::static_vector<size_t, 3>;

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <memory>
#include <vector>

struct TracklengthDistance {
  StaticVector3 index;
  double distance;
};

class PositionFilter {
 public:
  PositionFilter() = default;
  PositionFilter(std::size_t id) : id_(id) {}

  virtual ~PositionFilter() = default;

  virtual StaticVector3 get_indices(const Tracker& tktr) const = 0;

  virtual std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) const = 0;

  virtual StaticVector3 get_shape() const = 0;

  // method to get the volume and inverse of the volume
  virtual double dV(StaticVector3 index) const = 0;
  virtual double inv_dV(StaticVector3 index) const = 0;

  virtual std::string type_str() const = 0;

  std::size_t id() const { return id_; }

  virtual void write_to_hdf5(H5::Group& grp) const = 0;

 private:
  std::size_t id_;
};

std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node);

#endif