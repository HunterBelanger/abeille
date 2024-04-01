#ifndef GENERAL_TALLY_H
#define GENERAL_TALLY_H

#include <tallies/ITally.h>
#include <tallies/position_filter.h>
#include <tallies/energy_filter.h>

#include <tallies/box_position_filter.h>
#include <utils/mpi.hpp>


class GeneralTally : public ITally{
    public:
    GeneralTally(Quantity quantity, Estimator estimator, 
                std::shared_ptr<PositionFilter> position_filter, 
                std::shared_ptr<EnergyFilter> energy_in, std::shared_ptr<EnergyFilter> energy_out = nullptr)
        :ITally(quantity, estimator), position_filter_(position_filter),
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

    void score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat )override{
        std::array<int, 3> index_position;
        std::size_t index_E;

        if ( position_filter_->get_indices(tktr, index_position) 
        && energy_in_->get_index(p.E(), index_E) ){

            const double Et = mat.Et(p.E());
            double collision_score = 1.0 / ( Et * 1000);

            switch(quantity_){
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

        #ifdef ABEILLE_USE_OMP
        #pragma omp atomic
        #endif
            tally_gen_score( index_E,
                index_position[0], index_position[1], index_position[2]) += collision_score;
           
        }
    }

    private:
        std::shared_ptr<PositionFilter> position_filter_;
        std::shared_ptr<EnergyFilter> energy_in_;
        std::shared_ptr<EnergyFilter> energy_out_;
};

static Position Low_pos(-0.77032, -100., -100.);
static Position High_pos(0.77032, 100., 100.);

static std::shared_ptr<PositionFilter> const_pos = 
        std::make_shared<BoxPositionFilter>(Low_pos, High_pos);

static std::vector<double> energy_f{0.0, 1.0, 2.0};
static std::shared_ptr<EnergyFilter> const_ef
        = std::make_shared<EnergyFilter>(energy_f);

static std::shared_ptr<ITally> temp_tally
        = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux,
            GeneralTally::Estimator::Collision,
            const_pos, const_ef);


#endif