#include <tallies/box_position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/constants.hpp>

CartesianFilter::CartesianFilter(Position r_low, Position r_high)
    : r_low_(r_low), r_high_(r_high) {
  if ((r_low_.x() >= r_high_.x()) || (r_low_.y() >= r_high_.y()) ||
      (r_low_.z() >= r_high_.z()))
    fatal_error(
        " Coordinates of \"low\" position are >= than \"high\" "
        "position.\n");
}

// make the cartesian filter
std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node) {
  std::shared_ptr<CartesianFilter> cartesian_filter_ = nullptr;
  if (!node["position-filter"] && !node["position-filter"].IsScalar()) {
    fatal_error("position-filter is not given.");
  }
  const std::string cartesian_filter_type =
      node["position-filter"].as<std::string>();

  if (cartesian_filter_type == "box") {
    cartesian_filter_ = make_box_position_filter(node);
  } else if (cartesian_filter_type == "regular-cartesian-mesh") {
    cartesian_filter_ = make_regular_cartesian_filter(node);
  } else {
    fatal_error(cartesian_filter_type + " is not a valid position-filter.");
  }

  return cartesian_filter_;
}