#include <tallies/general_tally.hpp>
#include <utils/error.hpp>

void GeneralTally::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) {
    
    if (estimator_ == GeneralTally::Estimator::Collision){
        //std::array<int, 3> index_position;
        
        std::size_t index_E;

        if (energy_in_->get_index(p.E(), index_E) ){
        
            std::vector<size_t> indexes_;
            indexes_.reserve(4);

            indexes_ = position_filter_->get_indices(tktr);
            
            indexes_.insert(indexes_.begin(), index_E);
            
            if (indexes_.empty()){
                return;
            }

            

            const double Et = mat.Et(p.E());
            double collision_score = 1.0 / ( Et * net_weight_);

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
                tally_gen_score(indexes_) += collision_score;
        }   
    }
}

std::shared_ptr<ITally> temp_tally;

void make_temp_tally() {
    static bool init = false;

    if (!init) {
Position Low_pos(-0.77032, -100., -100.);
Position High_pos(0.77032, 100., 100.);

std::shared_ptr<PositionFilter> const_pos = 
    std::make_shared<BoxPositionFilter>(Low_pos, High_pos);

std::vector<double> energy_f{0.0, 1.0, 2.0};
std::shared_ptr<EnergyFilter> const_ef = std::make_shared<EnergyFilter>(energy_f);

temp_tally
    = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux,
        GeneralTally::Estimator::Collision, "temp_tally",
            const_pos, const_ef);
temp_tally->set_net_weight(1000.0);
      init = true;
    }

}



