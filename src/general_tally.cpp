#include <tallies/general_tally.h>
#include <utils/error.hpp>

void GeneralTally::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) {
    std::array<int, 3> index_position;
    std::size_t index_E;
    std::cout<<"Please come here.";
    fatal_error("Comes here.\n");

    if ( position_filter_->get_indices(tktr, index_position) 
    && energy_in_->get_index(p.E(), index_E) ){

        const double Et = mat.Et(p.E());
        double collision_score = 1.0 / ( Et * 1000);
        if (flag){
            fatal_error("Maybe coming.");
        }
        if (quantity_ == Quantity::Flux && flag3){
            fatal_error("BOOM.");
        }

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
            /*if (flag){
                std::cout<<"scored --->"<<tally_gen_score(index_E,
                index_position[0], index_position[1], index_position[2])<<"\n";
                flag = false;
            }  */
    }
}



