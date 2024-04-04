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
                size_t nx = position_filter_->Nx();
                size_t ny = position_filter_->Ny();
                size_t nz = position_filter_->Nz();
                size_t ne = energy_in_->size();
         
                tally_avg.reallocate({ne, nx, ny, nz});
                tally_avg.fill(0.0);

                tally_gen_score.reallocate({ne, nx, ny, nz});
                tally_gen_score.fill(0.0);

                tally_var.reallocate({ne, nx, ny, nz});
                tally_var.fill(0.0);
    }
    ~GeneralTally() = default;

    void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat )override;

    private:
        std::shared_ptr<PositionFilter> position_filter_;
        std::shared_ptr<EnergyFilter> energy_in_;
        std::shared_ptr<EnergyFilter> energy_out_;
};


extern std::shared_ptr<ITally> temp_tally;

#endif