#ifndef BOX_POSITION_FILTER_H
#define BOX_POSITION_FILTER_H

#include <array>

#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <yaml-cpp/yaml.h>

class BoxPositionFilter : public CartesianFilter{
    public:
    BoxPositionFilter(const Position r_low_, const Position r_high_)
    :CartesianFilter(r_low_, r_high_) {

    if ( (r_low.x() > r_high.x()) || (r_low.y() > r_high.y()) || (r_low.z() > r_high.z()) )
        fatal_error(" Corrdinates of \"low\" position are higher than \"high\" position.\n");
    
    dx = (r_high.x() - r_low.x());
    dy = (r_high.y() - r_low.y());
    dz = (r_high.z() - r_low.z());
    
    dx_inv = 1. / dx;
    dy_inv = 1. / dy;
    dz_inv = 1. / dz;

    }

    ~BoxPositionFilter() = default;

    /*bool get_indices(const Tracker& tktr, std::vector<std::size_t> indices)  {
        const Position r = tktr.r();
        if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
        && (r_low.y() <= r.y() && r_high.y() >= r.y())
        && (r_low.z() <= r.z() && r_high.z() >= r.z())){
            //indices.fill(static_cast<size_t> (0));
            indices.push_back(0);
            return true;
        }
        else{
            indices.push_back(-1);
            return false;
        }
    }*/

    std::vector<size_t> get_indices(const Tracker& tktr) override final{
        std::vector<size_t> indexes;
        indexes.reserve(1);
        const Position r = tktr.r();
        if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
        && (r_low.y() <= r.y() && r_high.y() >= r.y())
        && (r_low.z() <= r.z() && r_high.z() >= r.z())){
            
            indexes.push_back(0);

        }
        return indexes;
    }

    //std::vector<TracklengthDistance> get_indices_tracklength(Position r, const Direction& u_, double d_flight) override final;
    std::vector<TracklengthDistance> get_indices_tracklength(const Tracker& trkr, double d_flight) override final;

    void check_on_boundary(const Tracker tktr, std::array<int, 3> on);

    size_t Nx()const override final { return 1;}
    size_t Ny()const override final { return 1;}
    size_t Nz()const override final { return 1;}

    std::vector<size_t> get_dimension() override final{
        std::vector<size_t> pos_filter_dim{1};
        return pos_filter_dim;
    }

    double x_min(const std::vector<size_t>& index_)const override { 
        return r_low.x(); 
        if (index_.empty())
            return r_low.x();
    }
    double x_max(const std::vector<size_t>& index_)const override { 
        return r_high.x(); 
        if (index_.empty())
            return r_high.x();
    }
    
    double y_min(const std::vector<size_t>& index_)const override { 
        return r_low.y();
        if (index_.empty())
            return r_low.y();
    }
    double y_max(const std::vector<size_t>& index_)const override { 
        return r_high.y();
        if (index_.empty())
            return r_high.y();
    }
    
    double z_min(const std::vector<size_t>& index_)const override { 
        return r_low.z();
        if (index_.empty())
            return r_low.z();
    }

    double z_max(const std::vector<size_t>& index_)const override { 
        return r_high.z();
        if (index_.empty())
            return r_high.z();
    }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Box_Position_Filter; }
    std::string type_str()const override { return "box_position_filter"; }

};


//make the cartesian filter or position filter
template <typename BT>
std::shared_ptr<BT> make_box_position_filter(const YAML::Node &node){
    if (!node["low"])
        fatal_error("For box position-filter \"low\" co-ordinates is not provided.");
    if (!node["high"])
        fatal_error("For box position-filter \"high\" co-ordinates is not provided.");

    std::vector<double> low_point = node["low"].as<std::vector<double>>();
    std::vector<double> high_point = node["high"].as<std::vector<double>>();

    Position r_low_(low_point[0], low_point[1], low_point[2]);
    Position r_high_(high_point[0], high_point[1] , high_point[2]);

    std::shared_ptr<BT> box_type_filter 
                    = std::make_shared<BoxPositionFilter>(r_low_, r_high_);

    return box_type_filter;
}
#endif