#ifndef MESH_POSITION_FILTER_H
#define MESH_POSITION_FILTER_H


#include <tallies/cartesian_filter.hpp>
#include <utils/position.hpp>

#include <array>

class MeshPositionFilter : public CartesianFilter{
    public:
    MeshPositionFilter( Position r_low_, Position r_high_, 
        size_t nx_, size_t ny_, size_t nz_)
    : CartesianFilter(r_low_, r_high_),
    Nx_(nx_), Ny_(ny_), Nz_(nz_), 
    xmin(), ymin(), zmin()
    {
    
    if ( Nx_ == 0 || Ny_ == 0 || Nz_ == 0)
        fatal_error("The number of bins in any direction cannot be zero.\n");

    dx = (r_high.x() - r_low.x()) / static_cast<double>(Nx_);
    dy = (r_high.y() - r_low.y()) / static_cast<double>(Ny_);
    dz = (r_high.z() - r_low.z()) / static_cast<double>(Nz_);
    
    dx_inv = 1. / dx;
    dy_inv = 1. / dy;
    dz_inv = 1. / dz;
    

    }

    ~MeshPositionFilter() = default;

    std::vector<size_t> get_indices(const Tracker& tktr)override final;//override final;

    //std::vector<TracklengthDistance> get_indices_tracklength(Position r, const Direction& u_, double d_flight) override final;
    std::vector<TracklengthDistance> get_indices_tracklength(const Tracker& trkr, double d_flight) override final;
    
    size_t Nx()const override final { return Nx_;}
    size_t Ny()const override final { return Ny_;}
    size_t Nz()const override final { return Nz_;}

    double x_min()const override { return xmin; }
    double x_max()const override { return (xmin + dx); }

    double y_min()const override { return ymin; }
    double y_max()const override { return (ymin + dy); }
    
    double z_min()const override { return zmin; }
    double z_max()const override { return (zmin + dz); }

    std::vector<size_t> get_dimension() override final{
        std::vector<size_t> pos_filter_dim{Nx_, Ny_, Nz_};
        reduce_dimension(pos_filter_dim);
        return pos_filter_dim;
    }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Mesh_Positin_Filter; }
    std::string type_str()const override { return "Mesh_Position_Filter";}

    private:
    size_t Nx_, Ny_, Nz_;
    double xmin, ymin, zmin;

    //function will reduce the dimsion, if there is only one bin in the direction 
    void reduce_dimension(std::vector<size_t>& dimension_){
        int it_ = 0;
        if( (Nx_ == 1) && (dimension_.size() > 1) ){
            auto dimen_begin_ = dimension_.begin();
            dimension_.erase(dimen_begin_+ static_cast<std::ptrdiff_t>(it_));
            it_--;
        }

        it_++;
        if( (Ny_ == 1) && (dimension_.size() > 1) ){
            auto dimen_begin_ = dimension_.begin();
            dimension_.erase(dimen_begin_+ static_cast<std::ptrdiff_t>(it_));
            it_--;
        }

        it_++;
        if( (Nz_ == 1) && (dimension_.size() > 1) ){
            auto dimen_begin_ = dimension_.begin();
            dimension_.erase(dimen_begin_+ static_cast<std::ptrdiff_t>(it_));
            it_--;
        }

    }
};       

#endif