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

    std::vector<size_t> get_dimension() override final{
        std::vector<size_t> pos_filter_dim{Nx_, Ny_, Nz_};
        reduce_dimension(pos_filter_dim, 1);
        return pos_filter_dim;
    }

    //Perhaps Not Needed
    FilterType type()const override { return FilterType::Cartesian_Filter; }
    std::string type_str() const override { return "Cartesian_Filter"; };       

    protected:
        size_t Nx_, Ny_, Nz_;

        //function will reduce the dimsion, if there is only one bin in the direction 
        void reduce_dimension(std::vector<size_t>& dimension_, const std::size_t remove_num = 1){
            std::size_t it_ = 0;
            while(it_ < dimension_.size()){
                
                if( (dimension_[it_] == remove_num) && (dimension_.size() > 1) ){
                    auto dimen_begin_ = dimension_.begin();
                    dimension_.erase(dimen_begin_+ static_cast<std::ptrdiff_t>(it_));
                    continue;
                }
                it_++;
            }
        }
};

#endif