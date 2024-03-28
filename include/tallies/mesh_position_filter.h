#ifndef MESH_POSITION_FILTER_H
#define MESH_POSITION_FILTER_H

#include "position_filter.h"
#include <utils/position.hpp>

#include <array>

class Mesh_Position_Filter : public PositionFilter{
    public:
    Mesh_Position_Filter(Position r_low_, Position r_high_, int Nx_, int Ny_, int Nz_)
    : r_low(r_low_), r_high(r_high_), Nx(Nx_), Ny(Ny_), Nz(Nz_), 
    dx(), dy(), dz(), dx_inv(), dy_inv(), dz_inv() {
    
    if ( (r_low.x() > r_high.x()) && (r_low.y() > r_high.y()) && (r_low.z() > r_high.z()) )
        std::cout<<"Fatal Error, Corrdinates of \"low\" points are heigher than the corredinates of \"high\" points.\n";
    if ( Nx_ == 0 || Ny_ == 0 || Nz_ == 0)
        std::cout<<"Fatal Error, the number of bins in any direction should not be zero.\n";

    dx = (r_high.x() - r_low.x()) / static_cast<double>(Nx);
    dy = (r_high.y() - r_low.y()) / static_cast<double>(Ny);
    dz = (r_high.z() - r_low.z()) / static_cast<double>(Nz);
    
    dx_inv = 1. / dx;
    dy_inv = 1. / dy;
    dz_inv = 1. / dz;

    }

    ~Mesh_Position_Filter() = default;

    bool get_index(const Position& r, std::array<int, 3>& indices )override final {

        int index_x = static_cast<int> (std::floor( (r.x() - r_low.x()) * dx_inv) );
        int index_y = static_cast<int> (std::floor( (r.y() - r_low.y()) * dy_inv) );
        int index_z = static_cast<int> (std::floor( (r.z() - r_low.z()) * dz_inv) );

        if ( (index_x >= 0 && index_x < Nx) 
        && ( index_y >= 0 && index_y < Ny) 
        && (index_z >= 0 && index_z < Nz) ){
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
        
        return false;
    }

    FilterType type()const override { return FilterType::Mesh_Positin_Filter; }

    std::string type_str()const override { return "mesh_position_filter"; }

    int get_Nx()const { return Nx; }
    int get_Ny()const { return Ny; }
    int get_Nz()const { return Nz; }

    double x_min()const override { return xmin; }
    double x_max()const override { return (xmin + dx); }

    double y_min()const override { return ymin; }
    double y_max()const override { return (ymin + dy); }
    double z_min()const override { return zmin; }
    double z_max()const override { return (zmin + dz); }

    private:
    Position r_low, r_high;
    int Nx, Ny, Nz;
    double dx_inv, dy_inv, dz_inv, dx, dy, dz;
    double xmin, ymin, zmin;
};


#endif