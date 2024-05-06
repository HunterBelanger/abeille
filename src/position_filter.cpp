#include <tallies/box_position_filter.hpp>
#include <tallies/position_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <tallies/cylinder_position_filter.hpp>
#include <utils/error.hpp>

// make_position_filter will be usded for general tally system
std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node) {
  if (!node["position-filter"] && !node["position-filter"].IsScalar()) {
    fatal_error("position-filter is not given.");
  }

  const std::string position_filter_type =
      node["position-filter"].as<std::string>();

  std::shared_ptr<PositionFilter> position_filter_ = nullptr;

  if (position_filter_type == "box") {
    position_filter_ = make_box_position_filter(node);
  } else if (position_filter_type == "regular-cartesian-mesh") {
    position_filter_ = make_regular_cartesian_mesh_filter(node);
  } else if (position_filter_type == "cylinder-filter") {
    position_filter_ = make_cylinder_position_filter(node);
  } else {
    fatal_error(position_filter_type + " is not a valid position-filter.");
  }

  return position_filter_;
}