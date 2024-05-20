#include <tallies/cartesian_filter.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <utils/constants.hpp>

#include <sstream>

CartesianFilter::CartesianFilter(Position r_low, Position r_high,
                                 std::size_t id)
    : PositionFilter(id), r_low_(r_low), r_high_(r_high) {
  if ((r_low_.x() >= r_high_.x()) || (r_low_.y() >= r_high_.y()) ||
      (r_low_.z() >= r_high_.z())) {
    std::stringstream mssg;
    mssg << "CartesianFilter with id " << id
         << " does not satisfy r_low < r_high.";
    fatal_error(mssg.str());
  }
}

// make the cartesian filter
std::shared_ptr<CartesianFilter> make_cartesian_filter(const YAML::Node& node) {
  if (!node["type"] || node["type"].IsScalar() == false) {
    fatal_error("No valid type entry on position-filter.");
  }
  const std::string type = node["type"].as<std::string>();

  std::shared_ptr<CartesianFilter> cartesian_filter = nullptr;
  if (type == "regular-cartesian-mesh") {
    cartesian_filter = make_regular_cartesian_mesh_filter(node);
  } else {
    fatal_error("Unkown cartesian filter type " + type + ".");
  }

  return cartesian_filter;
}