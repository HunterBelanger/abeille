#ifndef CARTESIAN_FILTER_H
#define CARTESIAN_FILTER_H

#include <tallies/position_filter.hpp>
#include <utils/position.hpp>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>


class CartesianFilter : public PositionFilter{
    public:
    CartesianFilter(size_t nx_, size_t ny_, size_t nz_)
    : Nx_(nx_), Ny_(ny_), Nz_(nz_)
    {
    if ( Nx_ == 0 || Ny_ == 0 || Nz_ == 0)
        fatal_error("The number of bins in any direction cannot be zero.\n");
    }

    virtual ~CartesianFilter() = default;

    virtual double x_min()const = 0;
    virtual double x_max()const = 0;
    
    virtual double y_min()const = 0;
    virtual double y_max()const = 0;

    virtual double z_min()const = 0;
    virtual double z_max()const = 0;

    size_t Nx()const override { return Nx_; }
    size_t Ny()const override { return Ny_; }
    size_t Nz()const override { return Nz_; }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Cartesian_Filter; }
    std::string type_str() const override { return "Cartesian_Filter"; };       

    protected:
        size_t Nx_, Ny_, Nz_;
};

#endif