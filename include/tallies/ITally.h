#ifndef ITALLY_H
#define ITALLY_H

#include <memory>

#include "position_filter.h"
#include "energy_filter.h"
#include "ndarray.hpp"


class ITally{
    public: 
    enum class Quantity{
        Flux,
        Fission,
        Absorption,
        Elastic, 
    };

    enum class Axis{
        X, Y, Z
    };


    ITally(size_t FET_order_, Quantity quantity_, std::vector<Axis> axes_,
        std::shared_ptr<PositionFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) 
        :
        FET_coeff(), FET_var(), score_gen(), quantity(quantity_),
        axes(axes_), position_filter(position_filter_), 
        energy_in(energy_in_), energy_out(energy_out_)
        {
            //for legendre, check the size of axis vector should be between than 3.
            if ( axes.size() >= 1 && axes.size()<=3 )
                std::cout<<"Fatal error, the no. of given axis is not between 1 to 3.\n";
            else{
                
            }
            
            uint64_t Nx = static_cast<uint64_t> (position_filter->get_Nx());
            uint64_t Ny = static_cast<uint64_t> (position_filter->get_Ny());
            uint64_t Nz = static_cast<uint64_t> (position_filter->get_Nz());
            size_t ND  = axes.size();
            uint64_t NE = static_cast<uint64_t> (energy_in->get_size());
            // Order+1 is required as these many element is required. 
            size_t N_FET_element = static_cast<size_t> (FET_order_ + 1);

            FET_coeff.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_coeff.fill(0.0);

            score_gen.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            score_gen.fill(0.0);

            FET_var.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_var.fill(0.0);
            //FET_std.reallocate({NE, Nx, Ny, Nz});
            //FET_std.fill( UpperTriangularMatrix(FET_order + 1) );

            
        }
    
    ~ITally() = default;




     protected:
        std::shared_ptr<PositionFilter> position_filter {nullptr};
        std::shared_ptr<EnergyFilter> energy_in {nullptr};
        std::shared_ptr<EnergyFilter> energy_out {nullptr};
        
        NDArray<double> FET_coeff; // N-energy, N-x, N-y, N-z, N-axis for legendre, N_FET_order
        NDArray<double> score_gen;
        NDArray<double> FET_var;

        size_t gen = 0;

        Quantity quantity;
        std::vector<Axis> axes;
        
};


#endif