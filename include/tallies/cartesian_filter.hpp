#ifndef CARTESIAN_FILTER_H
#define CARTESIAN_FILTER_H

#include <simulation/tracker.hpp>
#include <tallies/position_filter.hpp>
#include <utils/error.hpp>
#include <utils/position.hpp>


class CartesianFilter : public PositionFilter {
 public:
  CartesianFilter(Position r_low_, Position r_high_)
      : r_low(r_low_), r_high(r_high_) {
    if ((r_low.x() > r_high.x()) || (r_low.y() > r_high.y()) ||
        (r_low.z() > r_high.z()))
      fatal_error(
          " Corrdinates of \"low\" position are higher than \"high\" "
          "position.\n");
  }

  virtual ~CartesianFilter() = default;

  virtual double x_min(const StaticVector6& index_) const = 0;
  virtual double x_max(const StaticVector6& index_) const = 0;

  virtual double y_min(const StaticVector6& index_) const = 0;
  virtual double y_max(const StaticVector6& index_) const = 0;

  virtual double z_min(const StaticVector6& index_) const = 0;
  virtual double z_max(const StaticVector6& index_) const = 0;

  virtual size_t Nx() const = 0;
  virtual size_t Ny() const = 0;
  virtual size_t Nz() const = 0;

  // required for track-length
  bool find_entry_point(Position& r, const Direction& u,
                        double& d_flight) const;
  void initialize_indices(const Position& r, const Direction& u, int& i, int& j,
                          int& k, std::array<int, 3>& on);
  void update_indices(int key, int& i, int& j, int& k, std::array<int, 3>& on);
  std::pair<double, int> distance_to_next_index(const Position& r,
                                                const Direction& u,
                                                const std::array<int, 3>& on,
                                                int i, int j, int k);

  // Perhaps Not Needed
  FilterType type() const override { return FilterType::Cartesian_Filter; }
  std::string type_str() const override { return "Cartesian_Filter"; };

 protected:
  Position r_low, r_high;
  double dx_inv, dy_inv, dz_inv, dx, dy, dz;
};

std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node);

#endif