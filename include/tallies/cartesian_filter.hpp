#ifndef CARTESIAN_FILTER_H
#define CARTESIAN_FILTER_H

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>

class CartesianFilter : public PositionFilter {
 public:
  CartesianFilter(Position r_low, Position r_high, std::size_t id);

  virtual ~CartesianFilter() = default;

  virtual StaticVector3 get_position_index(const Position& r) const = 0;

  virtual double x_min(const StaticVector3& index) const = 0;
  virtual double x_max(const StaticVector3& index) const = 0;
  virtual double dx(const StaticVector3& index) const = 0;
  virtual double inv_dx(const StaticVector3& index) const = 0;

  virtual double y_min(const StaticVector3& index) const = 0;
  virtual double y_max(const StaticVector3& index) const = 0;
  virtual double dy(const StaticVector3& index) const = 0;
  virtual double inv_dy(const StaticVector3& index) const = 0;

  virtual double z_min(const StaticVector3& index) const = 0;
  virtual double z_max(const StaticVector3& index) const = 0;
  virtual double dz(const StaticVector3& index) const = 0;
  virtual double inv_dz(const StaticVector3& index) const = 0;

  std::string type_str() const override { return "Cartesian_Filter"; };

 protected:
  Position r_low_, r_high_;
};

std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node);

#endif