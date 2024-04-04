#ifndef BOX_POSITION_FILTER_H
#define BOX_POSITION_FILTER_H

#include <array>

#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>


class BoxPositionFilter : public CartesianFilter{
    public:
    BoxPositionFilter(const Position r_low_, const Position r_high_)
    :CartesianFilter(1, 1, 1), r_low(r_low_), r_high(r_high_) {

    if ( (r_low.x() > r_high.x()) || (r_low.y() > r_high.y()) || (r_low.z() > r_high.z()) )
        fatal_error(" Corrdinates of \"low\" position are higher than \"high\" position.\n");
    }

    ~BoxPositionFilter() = default;

    bool get_indices(const Tracker& tktr, std::array<int, 3>& indices)  {
        const Position r = tktr.r();
        if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
        && (r_low.y() <= r.y() && r_high.y() >= r.y())
        && (r_low.z() <= r.z() && r_high.z() >= r.z())){
            indices.fill(static_cast<size_t> (0));
            return true;
        }
        else{
            indices.fill(-1);
            return false;
        }
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

    private:
    Position r_low, r_high;
};


#endif