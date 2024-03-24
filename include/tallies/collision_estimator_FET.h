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
    enum PolyEvalPlane{
        X, Y, Z
    };

    CollisionEstimator_FET(int FET_order_, Quantity quanitity_,
        std::shared_ptr<PositionFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) 
        :
        FET_order(FET_order_), FET_coeff(), FET_std(), quanitity(quanitity_),
        position_filter(position_filter_), 
        energy_in(energy_in_), energy_out(energy_out_)
        {
            uint64_t Nx = static_cast<uint64_t> (position_filter->get_Nx());
            uint64_t Ny = static_cast<uint64_t> (position_filter->get_Ny());
            uint64_t Nz = static_cast<uint64_t> (position_filter->get_Nz());
            uint64_t NE = static_cast<uint64_t> (energy_in->get_size()); 

            FET_coeff.reallocate({NE, Nx, Ny, Nz});
            FET_coeff.fill( std::vector<double> (FET_order_ + 1, 0.0));// Order+1 is required as these many element is required. 

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
            std::vector<double> temp_p_scr;
            temp_p_scr.resize( FET_order + 1);
            double beta_n;
            double x, y, z;
            for (int i = 0; i<FET_order+1; i++){
                
                switch(poly_type){
                    case PolynomialType::legendre:
                        x = 2 * ( x - position_filter->x_min()) / ( position_filter->x_min() - position_filter->x_max()) - 1 ;
                        beta_n = collision_score * legendre(i, x);
                        FET_coeff( index_E,
                        index_position[0], index_position[1], index_position[2])[i] += beta_n;
                        break;
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
        
        NDArray<std::vector<double>> FET_coeff;
        NDArray<UpperTriangularMatrix> FET_std;

        std::vector<std::vector<double>> particle_score; // for calculating the std deviation after particle death
        std::vector<std::vector<int>> particle_bin_index; // for storing the index

        int FET_order;  // Note that the number of coefficient will be one more than the order.
        Quantity quanitity; 
        PolynomialType poly_type = PolynomialType::legendre;

};


#endif