#ifndef CARTESIAN_FILTER_H
#define CARTESIAN_FILTER_H

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>

class CartesianFilter : public PositionFilter {
 public:
  CartesianFilter(Position r_low, Position r_high);

  virtual ~CartesianFilter() = default;

  virtual double x_min(const StaticVector3& index) const = 0;
  virtual double x_max(const StaticVector3& index) const = 0;

  virtual double y_min(const StaticVector3& index) const = 0;
  virtual double y_max(const StaticVector3& index) const = 0;

  virtual double z_min(const StaticVector3& index) const = 0;
  virtual double z_max(const StaticVector3& index) const = 0;

  virtual size_t Nx() const = 0;
  virtual size_t Ny() const = 0;
  virtual size_t Nz() const = 0;

  std::string type_str() const override { return "Cartesian_Filter"; };

 protected:

  Position r_low_, r_high_;

};

std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node);

#endif