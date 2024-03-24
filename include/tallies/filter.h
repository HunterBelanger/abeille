#ifndef FILTER_H
#define FILTER_H

#include <iostream>

//FilterType
enum class FilterType{
    Energy_Filter,
    Position_Filter,
    Box_Position_Filter,
    Mesh_Positin_Filter,
    Cylinder_Position_Filter,
    Cylinder_Array_Filter
};


//Filter class 
class Filter{
    public:
    Filter() = default;
    
    virtual ~Filter() = default;

    virtual FilterType type() const = 0;

    virtual std::string type_str()const = 0;
};

#endif

