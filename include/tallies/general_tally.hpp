#ifndef GENERAL_TALLY_H
#define GENERAL_TALLY_H

#include <tallies/itally.hpp>
#include <tallies/position_filter.hpp>
#include <tallies/energy_filter.hpp>

#include <tallies/box_position_filter.hpp>
#include <utils/mpi.hpp>

void make_temp_tally();

class GeneralTally : public ITally{
    public:
    GeneralTally(Quantity quantity, Estimator estimator, std::string name_,
                std::shared_ptr<PositionFilter> position_filter, 
                std::shared_ptr<EnergyFilter> energy_in, std::shared_ptr<EnergyFilter> energy_out = nullptr)
        :ITally(quantity, estimator, name_), position_filter_(position_filter),
        energy_in_(energy_in), energy_out_(energy_out)
            {
                /*size_t nx = position_filter_->Nx();
                size_t ny = position_filter_->Ny();
                size_t nz = position_filter_->Nz();
                
                tally_avg.reallocate({ne, nx, ny, nz});
                tally_avg.fill(0.0);

                tally_gen_score.reallocate({ne, nx, ny, nz});
                tally_gen_score.fill(0.0);

                tally_var.reallocate({ne, nx, ny, nz});
                tally_var.fill(0.0);*/

                //std::cout<<"---+++--------------------\n";
                const std::vector<size_t> dimen = position_filter_->get_dimension();
                
                //std::cout<<"\n-----------------------\n";
                // New constructor defination
                std::vector<size_t> tally_dimensions_;
                tally_dimensions_.reserve(4);
                tally_dimensions_ = position_filter_->get_dimension();
                size_t ne = energy_in_->size();
                tally_dimensions_.insert(tally_dimensions_.begin(), ne);

                
                const std::vector<size_t> temp_it = tally_dimensions_;
                tally_avg.reallocate(temp_it);
                tally_avg.fill(0.0);

                tally_gen_score.reallocate(temp_it);
                tally_gen_score.fill(0.0);

                tally_var.reallocate(temp_it);
                tally_var.fill(0.0);
                //std::cout<<"\n-----------------------\n";

    }
    ~GeneralTally() = default;

    void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) override final;

    //void score_flight(const Particle& p, double d_flight ,MaterialHelper& mat ) override final;
    void score_flight(const Particle& p, const Tracker& trkr, double d_flight, MaterialHelper& mat) override final;

    private:
        std::shared_ptr<PositionFilter> position_filter_;
        std::shared_ptr<EnergyFilter> energy_in_;
        std::shared_ptr<EnergyFilter> energy_out_;
};


extern std::shared_ptr<ITally> temp_tally;

#endif