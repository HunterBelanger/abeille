#include <tallies/cylinder_filter.hpp>
#include <tallies/position_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/error.hpp>

// make_position_filter will be usded for general tally system
std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node) {
  if (!node["type"] || node["type"].IsScalar() == false) {
    fatal_error("No type entry on position-filter.");
  }
  const std::string type = node["type"].as<std::string>();

  std::shared_ptr<PositionFilter> position_filter = nullptr;
  if (type == "regular-cartesian-mesh") {
    position_filter = make_regular_cartesian_mesh_filter(node);
  } else if (type == "cylinder-filter") {
    position_filter = make_cylinder_filter(node);
  } else {
    fatal_error("Unkown position-filter type " + type + ".");
  }

  return position_filter;
}