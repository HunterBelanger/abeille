#include <tallies/position_filter.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/box_position_filter.hpp>
#include <tallies/mesh_position_filter.hpp>


#include <utils/error.hpp>

// make_position_filter will be usded for general tally system
std::shared_ptr<PositionFilter> make_position_filter(const YAML::Node& node){
    std::shared_ptr<PositionFilter> position_filter_;

    if (!node["Position-Filter"]){
        fatal_error("Position-Filter is not given.");
    }
    const std::string position_filter_type = node["Position-Filter"].as<std::string>();

    if ( position_filter_type == "box"){
        position_filter_ = make_box_position_filter(node);
    }

    return position_filter_;
}