#include <tallies/box_position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/constants.hpp>

CartesianFilter::CartesianFilter(Position r_low, Position r_high)
    : r_low_(r_low), r_high_(r_high) {
  if ((r_low.x() > r_high.x()) || (r_low.y() > r_high.y()) ||
      (r_low.z() > r_high.z()))
    fatal_error(
        " Corrdinates of \"low\" position are higher than \"high\" "
        "position.\n");
}

// make the cartesian filter
std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node) {
  std::shared_ptr<CartesianFilter> cartesian_filter_ = nullptr;

  if (!node["Position-Filter"]) {
    fatal_error("Position-Filter is not given.");
  }
  const std::string cartesian_filter_type =
      node["Position-Filter"].as<std::string>();

  if (cartesian_filter_type == "box") {
    cartesian_filter_ = make_box_position_filter(node);
  }

  if (cartesian_filter_type == "regular-reactangular-mesh") {
    cartesian_filter_ = make_mesh_position_filter(node);
  }

  return cartesian_filter_;
}