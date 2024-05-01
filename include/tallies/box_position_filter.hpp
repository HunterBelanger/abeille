#ifndef BOX_POSITION_FILTER_H
#define BOX_POSITION_FILTER_H

#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

#include <array>

class BoxPositionFilter : public CartesianFilter {
 public:
  BoxPositionFilter(const Position r_low, const Position r_high);

  ~BoxPositionFilter() = default;

  StaticVector3 get_indices(const Tracker& tktr) override final {
    StaticVector3 indexes;
    const Position r = tktr.r();
    if ((r_low_.x() <= r.x() && r_high_.x() >= r.x()) &&
        (r_low_.y() <= r.y() && r_high_.y() >= r.y()) &&
        (r_low_.z() <= r.z() && r_high_.z() >= r.z())) {
      indexes.push_back(0);
    }
    return indexes;
  }

  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) override final;

  void check_on_boundary(const Tracker tktr, std::array<int, 3> on);

  size_t Nx() const override final { return 1; }
  size_t Ny() const override final { return 1; }
  size_t Nz() const override final { return 1; }

  StaticVector3 get_dimension() override final {
    return {1};
  }

  double x_min(const StaticVector3& /*index*/) const override {
    return r_low_.x();
  }
  double x_max(const StaticVector3& /*index*/) const override {
    return r_high_.x();
  }

  double y_min(const StaticVector3& /*index*/) const override {
    return r_low_.y();
  }
  double y_max(const StaticVector3& /*index*/) const override {
    return r_high_.y();
  }

  double z_min(const StaticVector3& /*index*/) const override {
    return r_low_.z();
  }

  double z_max(const StaticVector3& /*index*/) const override {
    return r_high_.z();
  }
  
  std::string type_str() const override { return "box_position_filter"; }

  protected:

    // required for track-length
  bool find_entry_point(Position& r, const Direction& u,
                        double& d_flight) const;
  void initialize_indices(const Position& r, const Direction& u, int& i, int& j,
                          int& k, std::array<int, 3>& on);
  std::pair<double, int> distance_to_next_index(const Position& r,
                                                const Direction& u,
                                                const std::array<int, 3>& on,
                                                int i, int j, int k);

  private:
    double dx_, dy_, dz_, dx_inv_, dy_inv_, dz_inv_ ;
};

// make the cartesian filter or position filter
std::shared_ptr<BoxPositionFilter> make_box_position_filter(const YAML::Node& node);

#endif