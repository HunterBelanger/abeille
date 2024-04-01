#ifndef ITALLY_H
#define ITALLY_H

#include <memory>

#include <ndarray.hpp>

#include <tallies/position_filter.h>
#include <tallies/energy_filter.h>
#include <utils/mpi.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <materials/material_helper.hpp>
#include <utils/error.hpp>



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

    enum class Estimator{
        Collision,
        TrackLength
    };

    ITally(Quantity quantity, Estimator estimator)
        :   quantity_(quantity), estimator_(estimator)
        {}
    
    virtual ~ITally() = default;

    virtual void score_tally(const Particle& p, const Tracker& tktr, MaterialHelper& mat){
        if(start_scoring_ == true){

            if ( estimator_ == Estimator::Collision){
                score_collision(p, tktr, mat );
            }
            //if ( estimator_ == Estimator::TrackLength){}
        }
    }

    //For collision estmator
    virtual void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) = 0;

    //For track-length estimator
    //virtual score_flight(const Particle& p, const Tracker& tktr, MaterialHelper& mat );

    //Record the avg and variance for the generation
    void record_tally(double mulitplier = 1.0);

    void set_net_weight(double weight_){ net_weight_ = weight_;}

    void start_scoring() { start_scoring_ = true;}
     //void end_scoring() { start_scoring_ = false;}

    protected:
        
        NDArray<double> tally_avg;
        NDArray<double> tally_gen_score;
        NDArray<double> tally_var;

        size_t gen_ = 0;
        double net_weight_;

        Quantity quantity_;
        Estimator estimator_;
    private:
        bool start_scoring_ = false;
        
        
};


#endif