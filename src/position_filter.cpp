#include <tallies/box_position_filter.hpp>
#include <tallies/position_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/error.hpp>

// make_position_filter will be usded for general tally system
std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node) {
  std::shared_ptr<PositionFilter> position_filter_ = nullptr;

  if (!node["Position-Filter"]) {
    std::string name_ = node["name"].as<std::string>();
    warning("Position-Filter is not given for the \"" + name_ + "\".");
    return position_filter_;
  }
  const std::string position_filter_type =
      node["Position-Filter"].as<std::string>();

  if (position_filter_type == "box") {
    position_filter_ = make_box_position_filter(node);
  }

  if (position_filter_type == "regular-reactangular-mesh") {
    position_filter_ = make_mesh_position_filter(node);
  }

  return position_filter_;
}