#ifndef FUNCTIONAL_EXPANSION_TALLIES_H
#define FUNCTIONAL_EXPANSION_TALLIES_H

#include <iostream>
#include <vector>
#include <memory>

#include "/Users/parthsingh/MonteCarlo/abeille/build/_deps/ndarray-src/include/ndarray.hpp"
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
        flux,
        fission,
        absorption,
        elastic, 
    };

    enum class PolynomialType{
        legendre,
        zernike
    };
    enum class Axis{
        X, Y, Z
    };

    CollisionEstimator_FET(int FET_order_, Quantity quanitity_,
        std::vector<Axis> axes_,
        std::shared_ptr<PositionFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) 
        :
        FET_order(FET_order_), FET_coeff(), FET_std(), quanitity(quanitity_),
        axes(axes_), position_filter(position_filter_), 
        energy_in(energy_in_), energy_out(energy_out_)
        {
            //for legendre, check the size of axis vector should be between than 3.
            if ( axes.size() >= 1 && axes.size()<=3 )
                std::cout<<"Fatal error, the no. of given axis is not between 1 to 3.\n";

            
            uint64_t Nx = static_cast<uint64_t> (position_filter->get_Nx());
            uint64_t Ny = static_cast<uint64_t> (position_filter->get_Ny());
            uint64_t Nz = static_cast<uint64_t> (position_filter->get_Nz());
            size_t ND  = axes.size();
            uint64_t NE = static_cast<uint64_t> (energy_in->get_size());
            // Order+1 is required as these many element is required. 
            size_t N_FET_element = static_cast<size_t> (FET_order + 1);

            FET_coeff.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_coeff.fill(0.0);

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
                case Quantity::flux:
                    collision_score *= p.wgt();
                    break;

                case Quantity::fission:
                    collision_score *= p.wgt() * mat.Ef(p.E());
                    break;

                case Quantity::absorption:
                    collision_score *= p.wgt() * mat.Ea(p.E());
                    break;

                case Quantity::elastic:
                    collision_score *= p.wgt() * mat.Eelastic(p.E());
                    break;
            }
            
            // adding the scores
            
            double beta_n;
            double x, y, z;
            for ( auto& c : axes){
                std::vector<double> temp_p_scr;
                temp_p_scr.reserve( FET_order + 1);
                for (int i = 0; i<FET_order+1; i++){

                    if (c == Axis::X){
                        x = 2 * (tktr.r().x() - position_filter->x_min()) 
                        / ( position_filter->x_max() - position_filter->x_min() ) - 1;
                        beta_n = collision_score * legendre(x, i) 
                            * lengendre_orthonormalization(i, position_filter->x_min(), position_filter->x_max());

                        temp_p_scr.   

                    }
                }
            }
        }
    }

    void score_into_tally(){

        
    }
    double polynomial_calc(Position R){
        
        return;
    }
    int get_FET_order() { return FET_order; } 

    private:
        std::shared_ptr<PositionFilter> position_filter {nullptr};
        std::shared_ptr<EnergyFilter> energy_in {nullptr};
        std::shared_ptr<EnergyFilter> energy_out {nullptr};
        
        NDArray<double> FET_coeff;
        NDArray<UpperTriangularMatrix> FET_std;
        std::vector<Axis> axes;

        std::vector<std::vector<double>> particle_score_x; // in-x,for calculating the std deviation after particle death
        std::vector<std::vector<double>> particle_score_y; // in-y, for calculating the std deviation after particle death
        std::vector<std::vector<double>> particle_score_z; // in-z, for calculating the std deviation after particle death

        std::vector<std::vector<int>> particle_bin_index; // for storing the index

        int FET_order;  // Note that the number of coefficient will be one more than the order.
        Quantity quanitity; 
        PolynomialType poly_type = PolynomialType::legendre;

};


#endif