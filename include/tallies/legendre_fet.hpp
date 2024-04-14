#ifndef LEGENDRE_FET_H
#define LEGENDRE_FET_H

#include <iostream>
#include <vector>
#include <memory>

#include <tallies/itally.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/legendre.hpp>

#include <ndarray.hpp>

#include <utils/mpi.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <materials/material_helper.hpp>

class LegendreFET : public ITally{
    public: 
    
    enum class Axis{
        X, Y, Z
    };


    LegendreFET(size_t FET_order_, Quantity quantity_, Estimator estimator_, std::string name_,
        std::vector<Axis> axes_,
        std::shared_ptr<CartesianFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in,
        std::shared_ptr<EnergyFilter> energy_out = nullptr) ;

    

    ~LegendreFET() = default;

    void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat );

    void score_flight(const Particle& p, const Tracker& trkr, double d_flight, MaterialHelper& mat){

    }

    size_t get_fet_order() { return FET_order; } 


    private:
        size_t FET_order;  // Note that the number of coefficient will be one more than the order.
        std::vector<Axis> axes;
        std::shared_ptr<CartesianFilter> cartesian_filter_;
        std::shared_ptr<EnergyFilter> energy_in_;
        std::shared_ptr<EnergyFilter> energy_out_;
        

};


#endif
