#ifndef POSITION_FILTER_H
#define POSITION_FILTER_H

#include "filter.h"
#include <utils/position.hpp>

class PositionFilter : public Filter{
    public:
    PositionFilter () = default;

    virtual ~PositionFilter () = default;

    
    virtual bool get_index(const Position& r, std::array<int, 3>& indices) = 0;
    virtual int get_Nx()const = 0;
    virtual int get_Ny()const = 0;
    virtual int get_Nz()const = 0;

    virtual double x_min()const = 0;
    virtual double x_max()const = 0;
    
    virtual double y_min()const = 0;
    virtual double y_max()const = 0;

    virtual double z_min()const = 0;
    virtual double z_max()const = 0;

    virtual FilterType type() const override { return FilterType::Position_Filter; };
    virtual std::string type_str() const override = 0;
       
};

#endif