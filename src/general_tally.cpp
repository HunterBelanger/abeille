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


//void GeneralTally::score_flight(const Particle& p,
//                                double d_flight ,MaterialHelper& mat ){

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr, 
                                double d_flight ,MaterialHelper& mat ){


    if (estimator_ == GeneralTally::Estimator::TrackLength){
        std::size_t index_E;

        if (energy_in_->get_index(p.E(), index_E) ){

            double flight_score = 1.0 / ( net_weight_);

            switch(quantity_){
                case Quantity::Flux:
                    flight_score *= p.wgt();
                    break;

                case Quantity::Fission:
                    flight_score *= p.wgt() * mat.Ef(p.E());
                    break;

                case Quantity::Absorption:
                    flight_score *= p.wgt() * mat.Ea(p.E());
                    break;

                case Quantity::Elastic:
                    flight_score *= p.wgt() * mat.Eelastic(p.E());
                    break;
            }

            std::vector<TracklengthDistance> pos_indexes_ = position_filter_->get_indices_tracklength(trkr, d_flight);

            for (size_t iter = 0; iter < pos_indexes_.size(); iter++){
                std::vector<size_t> all_indexes_ = pos_indexes_[iter].indexes_;
                const double d_ = pos_indexes_[iter].distance_in_bin;

                all_indexes_.insert(all_indexes_.begin(), index_E);

                #ifdef ABEILLE_USE_OMP
                #pragma omp atomic
                #endif
                    tally_gen_score(all_indexes_) += d_ * flight_score;
            }

        }
    }
}


std::shared_ptr<ITally> make_general_tally(const YAML::Node &node){

    
    if( !node["name"] ){
        fatal_error("No valid name is provided.");
    }
    std::string general_tally_name = node["name"].as<std::string>();

    if ( !node["estimator"]){
        fatal_error("No, estimator is given for \"" + general_tally_name
            + "\", therefore, collision estimator will be taken, by default.");
    }
    std::string estimator_name = node["estimator"].as<std::string>();

    
    std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;

    if (node["energy-bounds"]){
        std::vector<double> energy_bounds = node["energy-bounds"].as<std::vector<double>>();
        energy_filter_ = std::make_shared<EnergyFilter>(energy_bounds);
    }

    std::cout<<"energy-filter is made.\n";

    std::shared_ptr<PositionFilter> position_filter_ 
                                    = make_position_filter(node);
    
    std::cout<<"Position-filter is made."<<position_filter_->type_str()<<"\n";

    std::cout<<"general tally is made is made."<<estimator_name<<"\n";

    std::shared_ptr<ITally> itally_genral_tally = nullptr;

    if ( estimator_name == "collision" ){
        std::cout<<"general tally is made is made."<<estimator_name<<"\n";
        itally_genral_tally = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux, 
                                        GeneralTally::Estimator::Collision, general_tally_name,
                                             position_filter_, energy_filter_);
        std::cout<<"general tally is made is made.\n";
    }

    if ( estimator_name == "track-length" )
        itally_genral_tally = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux, 
                                        GeneralTally::Estimator::TrackLength, general_tally_name,
                                             position_filter_, energy_filter_);

    if ( itally_genral_tally == nullptr )
        fatal_error("Incorrect \"estimator\" is given.");

    return itally_genral_tally;

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





