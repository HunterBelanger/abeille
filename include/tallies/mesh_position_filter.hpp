#ifndef MESH_POSITION_FILTER_H
#define MESH_POSITION_FILTER_H


#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <array>

class MeshPositionFilter : public CartesianFilter{
    public:
    MeshPositionFilter( Position r_low_, Position r_high_, 
        size_t nx_, size_t ny_, size_t nz_)
    : CartesianFilter(nx_, ny_, nz_),
    r_low(r_low_), r_high(r_high_), 
    dx_inv(), dy_inv(), dz_inv(),
    xmin(), ymin(), zmin()
    {

    if ( (r_low.x() > r_high.x()) || (r_low.y() > r_high.y()) || (r_low.z() > r_high.z()) )
        fatal_error(" Corrdinates of \"low\" position are higher than \"high\" position.\n");

    dx = (r_high.x() - r_low.x()) / static_cast<double>(Nx_);
    dy = (r_high.y() - r_low.y()) / static_cast<double>(Ny_);
    dz = (r_high.z() - r_low.z()) / static_cast<double>(Nz_);
    
    dx_inv = 1. / dx;
    dy_inv = 1. / dy;
    dz_inv = 1. / dz;
    

    }

    ~MeshPositionFilter() = default;

    bool get_indices(const Tracker& tktr, std::array<int, 3>& indices )override final;

    double x_min()const override { return xmin; }
    double x_max()const override { return (xmin + dx); }

    double y_min()const override { return ymin; }
    double y_max()const override { return (ymin + dy); }
    
    double z_min()const override { return zmin; }
    double z_max()const override { return (zmin + dz); }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Mesh_Positin_Filter; }
    std::string type_str()const override { return "Mesh_Position_Filter";}

    private:
    Position r_low, r_high;
    double dx_inv, dy_inv, dz_inv, dx, dy, dz;
    double xmin, ymin, zmin;
};


bool MeshPositionFilter::get_indices(const Tracker& tktr, 
                                    std::array<int, 3>& indices ){
    
    const Position r = tktr.r();
    int index_x = static_cast<int> (std::floor( (r.x() - r_low.x()) * dx_inv) );
    int index_y = static_cast<int> (std::floor( (r.y() - r_low.y()) * dy_inv) );
    int index_z = static_cast<int> (std::floor( (r.z() - r_low.z()) * dz_inv) );

    if ( (index_x >= 0 && index_x < static_cast<int>(Nx_) ) 
    && ( index_y >= 0 && index_y < static_cast<int>(Ny_) ) 
    && (index_z >= 0 && index_z < static_cast<int>(Nz_) )){
        indices[0] = index_x;
        indices[1] = index_y;
        indices[2] = index_z;

        xmin = index_x * dx;
        ymin = index_y * dy;
        zmin = index_z * dz;

        return true;
    }
    else{
        indices.fill(-1);
        return false;
    }
}        
        

#endif