#ifndef POSITION_FILTER_H
#define POSITION_FILTER_H

#include <utils/error.hpp>
#include <utils/position.hpp>
#include <simulation/tracker.hpp>

enum class FilterType{
    Energy_Filter,
    Position_Filter,
    Cartesian_Filter,
    Box_Position_Filter,
    Mesh_Positin_Filter,
    Cylinder_Position_Filter,
    Cylinder_Array_Filter
};

class PositionFilter{
    public:
    PositionFilter() = default;
    
    virtual ~PositionFilter () = default;

    virtual bool get_indices(const Tracker& tktr, std::array<int, 3>& indices) = 0;
    
    virtual size_t Nx()const = 0;
    virtual size_t Ny()const = 0;
    virtual size_t Nz()const = 0;
    
    //Perhaps Not Needed
    virtual FilterType type()const { return FilterType::Position_Filter; }
    virtual std::string type_str() const = 0;
    

};

#endif