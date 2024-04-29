#include <tallies/itally.hpp>
#include <tallies/general_tally.hpp>
#include <utils/error.hpp>

#include <boost/container/static_vector.hpp>

template <typename std::size_t N>
using StaticVector = boost::container::static_vector<std::size_t, N>;


GeneralTally::GeneralTally(std::shared_ptr<PositionFilter> position_filter, 
                std::shared_ptr<EnergyFilter> energy_in,
                Quantity quantity, Estimator estimator, std::string name_)
    :
    ITally(quantity, estimator, name_),
    position_filter_(position_filter),
    energy_in_(energy_in)
{

    std::vector<size_t> tally_dimensions_;
    StaticVector<4> tally_dimensions_1;
    tally_dimensions_.reserve(4);

    bool check_ifany_filter = true;
    
    if ( position_filter_ ){ 
        tally_dimensions_ = position_filter_->get_dimension();
        check_ifany_filter = false;
    }
    
    if ( energy_in_ ){
        size_t ne = energy_in_->size();
        if ( ne > 1){
            tally_dimensions_.insert(tally_dimensions_.begin(), ne);
        }

        if ( (ne == 1) && (check_ifany_filter == true) ){
            tally_dimensions_.insert(tally_dimensions_.begin(), ne);
        }

        check_ifany_filter = check_ifany_filter && false;
    }

    if ( check_ifany_filter == true ){
        fatal_error("for the " + name_ + " tally, no filter is provided.");
    }

    tally_avg.reallocate(tally_dimensions_);
    tally_avg.fill(0.0);

    tally_gen_score.reallocate(tally_dimensions_);
    tally_gen_score.fill(0.0);

    tally_var.reallocate(tally_dimensions_);
    tally_var.fill(0.0);

}


void GeneralTally::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ) {
    
    if ( !(estimator_ == GeneralTally::Estimator::Collision) ){
        return;
    }

    std::size_t index_E;
    std::vector<size_t> indexes_;
    indexes_.reserve(4);

    bool check_filter = true;
    if ( position_filter_ ){
        indexes_ = position_filter_->get_indices(tktr); // it will provde the reduce dimensions
        check_filter = false;
    }

    if ( energy_in_ ){
        if (energy_in_->get_index(p.E(), index_E)){
            // Only add the energy_index if size of energy_bouds are more than one.
            if ( energy_in_->size() > 1 ){
                indexes_.insert(indexes_.begin(), index_E);
            }

            if ( energy_in_->size() == 1 && check_filter == true){
                indexes_.insert(indexes_.begin(), index_E);
            }

        }
        else{
            return;
        }
    }
    
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


//void GeneralTally::score_flight(const Particle& p,
//                                double d_flight ,MaterialHelper& mat ){

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr, 
                                double d_flight ,MaterialHelper& mat ){


    if ( !(estimator_ == GeneralTally::Estimator::TrackLength) ){
        return;
    }
        
    std::size_t index_E;
    bool energy_in_multi_size = false;  // Should be true when the energy_in_->size() is greater than 1.

    if ( energy_in_ ){
        if (energy_in_->get_index(p.E(), index_E) ){
            // Only add the energy_index if size of energy_bouds are more than one.
            // for track-length, it will done with help of bool
            if ( energy_in_->size() > 1)
                energy_in_multi_size = true;
        }
        else{
            return;
        }
    }

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

        if (energy_in_multi_size == true){
            all_indexes_.insert(all_indexes_.begin(), index_E);
        }

        #ifdef ABEILLE_USE_OMP
        #pragma omp atomic
        #endif
            tally_gen_score(all_indexes_) += d_ * flight_score;
    }
}


std::shared_ptr<ITally> make_general_tally(const YAML::Node &node){

    // Check the name of the tally is given or not.
    if( !node["name"] ){
        fatal_error("No valid name is provided.");
    }
    std::string general_tally_name = node["name"].as<std::string>();

    // Check the estimator is given or not.
    if ( !node["estimator"]){
        fatal_error("No, estimator is given for \"" + general_tally_name
            + "\", therefore, collision estimator will be taken, by default.");
    }
    std::string estimator_name = node["estimator"].as<std::string>();

    // Get the enrgy bounds, if any is given
    std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;
    if (node["energy-bounds"]){
        std::vector<double> energy_bounds = node["energy-bounds"].as<std::vector<double>>();
        energy_filter_ = std::make_shared<EnergyFilter>(energy_bounds);
    }

    // Get the position filter
    std::shared_ptr<PositionFilter> position_filter_ 
                                    = make_position_filter(node);

    // For the general tally
    std::shared_ptr<ITally> itally_genral_tally = nullptr;
    if ( estimator_name == "collision" ){
        itally_genral_tally = std::make_shared<GeneralTally>(position_filter_, energy_filter_,
                                                    GeneralTally::Quantity::Flux, 
                                        GeneralTally::Estimator::Collision, general_tally_name);
    }

    if ( estimator_name == "track-length" )
        itally_genral_tally = std::make_shared<GeneralTally>(position_filter_, energy_filter_,
                                                    GeneralTally::Quantity::Flux, 
                                        GeneralTally::Estimator::TrackLength, general_tally_name);

    if ( itally_genral_tally == nullptr )
        fatal_error("Incorrect \"estimator\" is given.");

    return itally_genral_tally;

}
