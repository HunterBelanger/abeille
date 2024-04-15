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

    double x_min()const override { return r_low.x(); }
    double x_max()const override { return r_high.x(); }
    
    double y_min()const override { return r_low.y(); }
    double y_max()const override { return r_high.y(); }
    
    double z_min()const override { return r_low.z(); }
    double z_max()const override { return r_high.z(); }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Box_Position_Filter; }
    std::string type_str()const override { return "box_position_filter"; }

};


//make the cartesian filter or position filter
//template<typename T>
std::shared_ptr<PositionFilter> make_box_position_filter(const YAML::Node &node);


#endif