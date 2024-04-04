#ifndef LEGENDRE_COLLISION_FET_H
#define LEGENDRE_COLLISION_FET_H

#include <iostream>
#include <vector>
#include <memory>

#include "ITally.h"

#include "ndarray.hpp"
#include "box_position_filter.h"
#include "mesh_position_filter.h"
#include "cylinder_position_filter.h"
#include "upper_triangular_matrix.h"
#include "energy_filter.h"
#include "legendre.h"

#include <utils/mpi.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <materials/material_helper.hpp>

class CollisionEstimatorFET : public ITally{
    public: 
    CollisionEstimator_FET(size_t FET_order_, Quantity quantity_,
        std::vector<Axis> axes_,
        std::shared_ptr<PositionFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) ;

    
    enum class Axis{
        X, Y, Z
    };
    ~CollisionEstimatorFET() = default;

    void score_particle(const Particle& p, const Tracker& tktr, MaterialHelper& mat );
    
    void record_variance(double multiplier);

    double polynomial_calc(Position R){
        
        return;
    }
    
    size_t get_FET_order() { return FET_order; } 


    private:
        size_t FET_order;  // Note that the number of coefficient will be one more than the order.
        std::vector<Axis> axes;
        

};


#endif







/*
            for (size_t it_E = 0; it_E < energy_in->get_size(); it_E++){
                for(size_t it_x = 0; it_x < position_filter->get_Nx(); it_x++){
                    for (size_t it_y = 0; it_y < position_filter->get_Ny(); it_y++){
                        for (size_t it_z = 0; it_z < position_filter->get_Nz(); it_z++){
                            for (size_t it_axis = 0; it_axis < axes.size(); it_axis++ ){
                                for (size_t it_coeff = 0; it_coeff < FET_order+1; it_coeff++){
                                    
                                    
                                    
                                }
                            }
                        }
                    }
                }
            }

*/