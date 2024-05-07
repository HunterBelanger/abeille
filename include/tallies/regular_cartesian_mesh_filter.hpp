#ifndef REGULAR_CARTESIAN_MESH_FILTER_H
#define REGULAR_CARTESIAN_MESH_FILTER_H

#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <array>

class RegularCartesianMeshFilter : public CartesianFilter {
 public:
  RegularCartesianMeshFilter(Position r_low, Position r_high, size_t nx,
                             size_t ny, size_t nz, std::size_t id);

  StaticVector3 get_indices(const Tracker& tktr) override final;

  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) override final;

  double x_min(const StaticVector3& index) const override;
  double x_max(const StaticVector3& index) const override;
  double y_min(const StaticVector3& index) const override;
  double y_max(const StaticVector3& index) const override;
  double z_min(const StaticVector3& index) const override;
  double z_max(const StaticVector3& index) const override;

  StaticVector3 get_shape() override final {
    if (Nx_ == 1 && Ny_ == 1 && Nz_ == 1){
      return {1};
    }
    return reduce_dimension(Nx_, Ny_, Nz_);
  }

  std::string type_str() const override { return "regular-cartesian-mesh"; }

 protected:
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

 private:
  double dx_, dy_, dz_, dx_inv_, dy_inv_, dz_inv_;
  size_t Nx_, Ny_, Nz_, x_index_, y_index_, z_index_;

  // function will reduce the dimsion, if there is only one bin in the direction
  StaticVector3 reduce_dimension(const size_t& loc_x, const size_t& loc_y,
                                 const size_t& loc_z) {
    if (Nx_ == 1 && Ny_ == 1 && Nz_ == 1){
      return {loc_x};
    }
    StaticVector3 reduce_;
    if (Nx_ > 1) {
      reduce_.push_back(loc_x);
    }

    if (Ny_ > 1) {
      reduce_.push_back(loc_y);
    }

    if (Ny_ > 1) {
      reduce_.push_back(loc_z);
    }
    return reduce_;
  }
};

// Make the cartesian or position filter class
std::shared_ptr<RegularCartesianMeshFilter> make_regular_cartesian_mesh_filter(
    const YAML::Node& node);

#endif