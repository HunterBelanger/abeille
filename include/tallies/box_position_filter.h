#ifndef BOX_POSITION_FILTER_H
#define BOX_POSITION_FILTER_H

#include <array>

#include "position_filter.h"
#include <utils/position.hpp>


class BoxPositionFilter : public PositionFilter{
    public:
    BoxPositionFilter(const Position& r_low_, const Position& r_high_)
    : r_low(r_low_), r_high(r_high_) {}

    ~BoxPositionFilter() = default;

    bool get_index(const Position& r, std::array<int, 3>& indices)override final {
        if ( (r_low.x() <= r.x() && r_high.x() >= r.x())
        && (r_low.y() <= r.y() && r_high.y() >= r.y())
        && (r_low.z() <= r.z() && r_high.z() >= r.z())){
            indices.fill(static_cast<int> (0));
            return true;
        }
        else{
            indices.fill(-1);
            return false;
        }
    }

    int get_Nx()const override { return 1; }
    int get_Ny()const override { return 1; }
    int get_Nz()const override { return 1; }

    double x_min()const override { return r_low.x(); }
    double x_max()const override { return r_high.x(); }
    double y_min()const override { return r_low.y(); }
    double y_max()const override { return r_high.y(); }
    double z_min()const override { return r_low.z(); }
    double z_max()const override { return r_high.z(); }

    FilterType type()const override { return FilterType::Box_Position_Filter; }

    std::string type_str()const override { return "box_position_filter"; }

    private:
    Position r_low, r_high;
};


#endif