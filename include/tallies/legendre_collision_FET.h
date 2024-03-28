#ifndef FUNCTIONAL_EXPANSION_TALLIES_H
#define FUNCTIONAL_EXPANSION_TALLIES_H

#include <iostream>
#include <vector>
#include <memory>

#include "ndarray.hpp"
#include "box_position_filter.h"
#include "mesh_position_filter.h"
#include "cylinder_position_filter.h"
#include "upper_triangular_matrix.h"
#include "energy_filter.h"
#include "legendre.h"

#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <materials/material_helper.hpp>

class CollisionEstimator_FET{
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

    CollisionEstimator_FET(size_t FET_order_, Quantity quanitity_,
        std::vector<Axis> axes_,
        std::shared_ptr<PositionFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) 
        :
        FET_order(FET_order_), FET_coeff(), FET_std(), score_gen(), quanitity(quanitity_),
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
            size_t N_FET_element = static_cast<size_t> (FET_order + 1);

            FET_coeff.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_coeff.fill(0.0);

            score_gen.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            score_gen.fill(0.0);

            FET_std.reallocate({NE, Nx, Ny, Nz});
            FET_std.fill( UpperTriangularMatrix(FET_order + 1) );

            
        }
    
    ~CollisionEstimator_FET() = default;

    void score_particle(const Particle& p, const Tracker& tktr, MaterialHelper& mat ){
        std::array<int, 3> index_position;
        std::size_t index_E;

        if ( position_filter->get_index(tktr.r(), index_position) 
        && energy_in->get_index(p.E(), index_E) ){

            const double Et = mat.Et(p.E());

            double collision_score = 1.0 / Et;

            switch(quanitity){
                case Quantity::Flux:
                    collision_score *= p.wgt();
                    break;

                case Quantity::Fission:
                    collision_score *= p.wgt() * mat.Ef(p.E());
                    break;

                case Quantity::Absorption:
                    collision_score *= p.wgt() * mat.Ea(p.E());
                    break;

                case Quantity::Elastic:
                    collision_score *= p.wgt() * mat.Eelastic(p.E());
                    break;
            }
            
            // adding the scores
            
            double beta_n;
            size_t it_axis = 0;
            for ( auto& c : axes){  // Loop over the different axis
                for (size_t i = 0; i<FET_order+1; i++){    // loop over differnt FET order
                    if (c == Axis::X){
                        const double xmin_ = position_filter->x_min();
                        const double xmax_ = position_filter->x_max();
                        const double x = 2 * (tktr.r().x() - xmin_) / (xmax_ - xmin_);

                        const double beta_n = collision_score * legendre(i ,x); // score for i-th order's basis function
                        FET_coeff(index_E, 
                                index_position[0], index_position[1], index_position[2],
                                it_axis, i );

                        
                        score_gen(index_E, 
                                index_position[0], index_position[1], index_position[2],
                                it_axis, i ) += beta_n;
                                

                        it_axis++;
                    }
                }
            }
        }
    }

    void record_variance(){
        gen++;
        const double  inv_gen = 1.0 / ( static_cast<double> (gen));

        for( auto &c : axes){
            for( size_t i = 0; i < FET_coeff.size(); i++){
                double old_avg = score_gen[i];
                double temp_value = FET_coeff[i];

                
            }
        }

        
    }

    double polynomial_calc(Position R){
        
        return;
    }
    
    size_t get_FET_order() { return FET_order; } 

    private:
        std::shared_ptr<PositionFilter> position_filter {nullptr};
        std::shared_ptr<EnergyFilter> energy_in {nullptr};
        std::shared_ptr<EnergyFilter> energy_out {nullptr};
        
        NDArray<double> FET_coeff; // N-energy, N-x, N-y, N-z, N-axis for legendre, N_FET_order
        NDArray<double> score_gen;
        NDArray<UpperTriangularMatrix> FET_std;
        
        std::vector<Axis> axes;

        size_t FET_order;  // Note that the number of coefficient will be one more than the order.
        Quantity quanitity; 

    protected:
        size_t gen = 0;
        

};


#endif