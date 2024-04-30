#ifndef MESH_POSITION_FILTER_H
#define MESH_POSITION_FILTER_H

#include <tallies/box_position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <array>

class MeshPositionFilter : public CartesianFilter {
 public:
  MeshPositionFilter(Position r_low_, Position r_high_, size_t nx_, size_t ny_,
                     size_t nz_)
      : CartesianFilter(r_low_, r_high_),
        Nx_(nx_),
        Ny_(ny_),
        Nz_(nz_),
        ymin(),
        zmin(),
        index_x(),
        index_y(),
        index_z() {
    if (Nx_ == 0 || Ny_ == 0 || Nz_ == 0)
      fatal_error("The number of bins in any direction cannot be zero.\n");

    dx = (r_high.x() - r_low.x()) / static_cast<double>(Nx_);
    dy = (r_high.y() - r_low.y()) / static_cast<double>(Ny_);
    dz = (r_high.z() - r_low.z()) / static_cast<double>(Nz_);

    dx_inv = 1. / dx;
    dy_inv = 1. / dy;
    dz_inv = 1. / dz;

    std::cout << "Dx = " << dx << "\n(1) " << r_low.x() << "\n(2)"
              << (r_low.x() + dx) << "\n(3)" << (r_low.x() + 2.0 * dx) << "\n";

    index_x = 0;
    index_y = 1;
    index_z = 2;
    if (Nx_ == 1){
        index_x = 0;
        index_y--;
        index_z--;
    }

    if ( Ny_ == 1){
        index_y = 0;
        index_z--;
  }

  if (Nz_ == 1){
    index_z = 0;
  }

  
        }

  ~MeshPositionFilter() = default;

  StaticVector6 get_indices(
      const Tracker& tktr) override final;  // override final;

  // std::vector<TracklengthDistance> get_indices_tracklength(Position r, const
  // Direction& u_, double d_flight) override final;
  std::vector<TracklengthDistance> get_indices_tracklength(
      const Tracker& trkr, double d_flight) override final;

  size_t Nx() const override final { return Nx_; }
  size_t Ny() const override final { return Ny_; }
  size_t Nz() const override final { return Nz_; }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double x_min(const StaticVector6& index_) const override {
    if (Nx_ == 1) {
      return r_low.x();
    }

    const double xmin_ = r_low.x() + static_cast<double>(index_[index_x]) * dx;
    return xmin_;
  }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double x_max(const StaticVector6& index_) const override {
    if (Nx_ == 1) {
      return r_high.x();
    }
    const double xmax_ = r_low.x() + static_cast<double>(index_[index_x]) * dx + dx;
    return xmax_;
  }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double y_min(const StaticVector6& index_) const override {
    if (Ny_ == 1) {
      return r_low.y();
    }

    const double ymin_ =
        r_low.y() + static_cast<double>(index_[index_y]) * dy;
    return ymin_;
  }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double y_max(const StaticVector6& index_) const override {
    if (Ny_ == 1) return r_high.y();

    const double ymax_ =
        r_low.y() + static_cast<double>(index_[index_y]) * dy + dy;
    return ymax_;
  }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double z_min(const StaticVector6& index_) const override {
    if (Nz_ == 1) return r_low.z();

    const double zmin_ =
        r_low.z() + static_cast<double>(index_[index_z]) * dz;
    return zmin_;
  }

  // Note that the method given is applicable when the reduced_dimension is
  // used.
  double z_max(const StaticVector6& index_) const override {
    if (Nz_ == 1) return r_high.z();

    const double zmax_ =
        r_low.z() + static_cast<double>(index_[index_z]) * dz + dz;
    return zmax_;
  }

  StaticVector6 get_dimension() override final {
    StaticVector6 pos_filter_dim{Nx_, Ny_, Nz_};
    reduce_dimension(pos_filter_dim);
    return pos_filter_dim;
  }

  // Perhaps Not Needed
  FilterType type() const override { return FilterType::Mesh_Positin_Filter; }
  std::string type_str() const override { return "Mesh_Position_Filter"; }

 private:
  size_t Nx_, Ny_, Nz_;
  double ymin, zmin;
  size_t index_x, index_y, index_z;

  // function will reduce the dimsion, if there is only one bin in the direction
  void reduce_dimension(StaticVector6& dimension_) {
    int it_ = 0;
    if ((Nx_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }

    it_++;
    if ((Ny_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }

    it_++;
    if ((Nz_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }
  }

    // function-overloading for std::vecotr<size_t> type will reduce the dimsion, if there is only one bin in the direction
    // this is currently being used for the track-length
  void reduce_dimension(std::vector<size_t> dimension_) {
    int it_ = 0;
    if ((Nx_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }

    it_++;
    if ((Ny_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }

    it_++;
    if ((Nz_ == 1) && (dimension_.size() > 1)) {
      auto dimen_begin_ = dimension_.begin();
      dimension_.erase(dimen_begin_ + static_cast<std::ptrdiff_t>(it_));
      it_--;
    }
  }


};

// Make the cartesian or position filter class
template <typename T>
std::shared_ptr<T> make_mesh_position_filter(const YAML::Node& node) {
  if (!node["low"])
    fatal_error(
        "For mesh position-filter \"low\" co-ordinates is not provided.");
  if (!node["high"])
    fatal_error(
        "For mesh position-filter \"high\" co-ordinates is not provided.");
  if (!node["shape"])
    fatal_error("For mesh position-filter \"shape\" is not provided.");

  std::vector<double> low_point = node["low"].as<std::vector<double>>();
  std::vector<double> high_point = node["high"].as<std::vector<double>>();
  std::vector<std::size_t> shape_ =
      node["shape"].as<std::vector<std::size_t>>();

  if (shape_[0] == 1 && shape_[1] == 1 && shape_[2] == 1) {
    std::string name_ = node["name"].as<std::string>();
    warning("For the tally " + name_ +
            ", given shape is [1, 1, 1] in reactangular-mesh, therefore, a "
            "\"box\" filter will be used.");
    std::shared_ptr<T> mesh_box_filter_ = make_box_position_filter<T>(node);
    return mesh_box_filter_;
  }

  if (shape_.size() != 3) fatal_error("Element in the shape must be 3.");

  Position r_low_(low_point[0], low_point[1], low_point[2]);
  Position r_high_(high_point[0], high_point[1], high_point[2]);

  std::shared_ptr<T> mesh_type_filter = std::make_shared<MeshPositionFilter>(
      r_low_, r_high_, shape_[0], shape_[1], shape_[2]);

  return mesh_type_filter;
}

#endif