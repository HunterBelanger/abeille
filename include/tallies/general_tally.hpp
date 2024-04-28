#ifndef GENERAL_TALLY_H
#define GENERAL_TALLY_H

#include <memory>

#include <yaml-cpp/yaml.h>

#include <tallies/itally.hpp>
#include <tallies/position_filter.hpp>
#include <tallies/energy_filter.hpp>

#include <tallies/box_position_filter.hpp>
#include <utils/mpi.hpp>



class GeneralTally : public ITally{
    public:
    GeneralTally(std::shared_ptr<PositionFilter> position_filter, 
                std::shared_ptr<EnergyFilter> energy_in,
                Quantity quantity, Estimator estimator, std::string name_);

    /*
    // The below defination can be used when energ_out is not a nullptr
    GeneralTally(std::shared_ptr<PositionFilter> position_filter, std::shared_ptr<EnergyFilter> energy_in, 
                    std::shared_ptr<EnergyFilter> energy_out, Quantity quantity, Estimator estimator, std::string name_);
    */

    ~GeneralTally() = default;

    void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) override final;

    //void score_flight(const Particle& p, double d_flight ,MaterialHelper& mat ) override final;
    void score_flight(const Particle& p, const Tracker& trkr, double d_flight, MaterialHelper& mat) override final;

    private:
        std::shared_ptr<PositionFilter> position_filter_;
        std::shared_ptr<EnergyFilter> energy_in_;
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr;


};

std::shared_ptr<ITally> make_general_tally(const YAML::Node &node);


#endif