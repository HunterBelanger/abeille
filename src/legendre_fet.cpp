#include <tallies/legendre_fet.hpp>
#include <utils/mpi.hpp>

#include <tallies/itally.hpp>

#include <boost/container/static_vector.hpp>



LegendreFET::LegendreFET(size_t FET_order_, LegendreFET::Quantity quantity_, LegendreFET::Estimator estimator_, std::string name_,
        std::vector<LegendreFET::Axis> axes_,
        std::shared_ptr<CartesianFilter> position_filter_, 
        std::shared_ptr<EnergyFilter> energy_in,
        std::shared_ptr<EnergyFilter> energy_out)
    :
    ITally(quantity_, estimator_, name_), 
    FET_order(FET_order_), axes(axes_), cartesian_filter_ (position_filter_), 
    energy_in_(energy_in), energy_out_(energy_out)
    {
    
    //for legendre, check the size of axis vector should be between than 3.
    if ( axes.size() <= 1 && axes.size() >=3 )
        fatal_error("The no. of given axis for legendre-FET is not between 1 to 3.");

    std::vector<size_t> tally_dimensions_ = cartesian_filter_->get_dimension();
    tally_dimensions_.reserve(6);
    std::size_t ne = energy_in_->size();
    tally_dimensions_.insert(tally_dimensions_.begin(), ne);
    tally_dimensions_.push_back(axes.size());
    tally_dimensions_.push_back(FET_order_+1);


    tally_avg.reallocate(tally_dimensions_);
    tally_avg.fill(0.0);

    tally_gen_score.reallocate(tally_dimensions_);
    tally_gen_score.fill(0.0);
    
    tally_var.reallocate(tally_dimensions_);
    tally_var.fill(0.0);


    std::cout<<"FET-Tally Dimemsion : "<<tally_dimensions_.size()<<"\n";
    

    }

void LegendreFET::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ){
    std::vector<size_t> index_position;
    std::size_t index_E;

    if ( energy_in_->get_index(p.E(), index_E) ){

        index_position = cartesian_filter_->get_indices(tktr) ;

        if (index_position.empty()){
            return;
        }

        const double Et = mat.Et(p.E());

        double collision_score = 1.0 / ( Et * net_weight_);

        switch(quantity_){
            case LegendreFET::Quantity::Flux:
                collision_score *= p.wgt();
                break;

            case LegendreFET::Quantity::Fission:
                collision_score *= p.wgt() * mat.Ef(p.E());
                break;

            case LegendreFET::Quantity::Absorption:
                collision_score *= p.wgt() * mat.Ea(p.E());
                break;

            case LegendreFET::Quantity::Elastic:
                collision_score *= p.wgt() * mat.Eelastic(p.E());
                break;
        }
        
        // adding the scores
        std::vector<size_t> indexes_;
        indexes_.reserve(6);
        indexes_ = cartesian_filter_->get_indices(tktr) ;

        if (index_position.empty()){
            return;
        }
            
        indexes_.insert(indexes_.begin(), index_E);
    

        const size_t axis_index = indexes_.size();
        const size_t FET_index = indexes_.size() + 1;
        
        indexes_.push_back(0);
        indexes_.push_back(0);

        double beta_n;

        boost::container::static_vector<std::size_t, 10> statc_vec;
        for ( auto &c : indexes_ ){
            statc_vec.push_back(c);
        }

        size_t it_axis = 0;
        //for ( auto& c : axes){  // Loop over the different axis and indexing is done using the it_axis
            for (size_t i = 0; i<FET_order+1; i++){    // loop over differnt FET order

               /* switch (c){
                case LegendreFET::Axis::X: {*/
                    const double xmin_ = cartesian_filter_->x_min(indexes_);
                    const double xmax_ = cartesian_filter_->x_max(indexes_);
                    const double x = 2. * (tktr.r().x() - xmin_) / (xmax_ - xmin_) - 1.;

                    beta_n = collision_score * legendre(i ,x); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
                    
            #ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                tally_gen_score(statc_vec) += beta_n;
            }
                    
            /*#ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;
            */
            /*       
                    }
                    break;
                
                case LegendreFET::Axis::Y: {
                    const double ymin_ = cartesian_filter_->y_min();
                    const double ymax_ = cartesian_filter_->y_max();
                    const double y = 2 * (tktr.r().y() - ymin_) / (ymax_ - ymin_) - 1;

                    beta_n = collision_score * legendre(i ,y); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
            /*#ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;    
            */
            /*        }
                    break;
                
                case LegendreFET::Axis::Z: {
                    const double zmin_ = cartesian_filter_->z_min();
                    const double zmax_ = cartesian_filter_->z_max();
                    const double z = 2 * (tktr.r().z() - zmin_) / (zmax_ - zmin_) - 1;

                    beta_n = collision_score * legendre(i , z); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
            /*#ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;

                    }
                    break;
                
                }*/



            it_axis++;
        //}

    }
}



std::shared_ptr<ITally> make_legendre_fet(const YAML::Node& node){

    // Check the name of the tally is given or not.
    if( !node["name"] ){
        fatal_error("No valid name is provided.");
    }
    std::string legendre_fet_tally_name = node["name"].as<std::string>();  

    // Check the estimator is given or not.
    if ( !node["estimator"]){
        fatal_error("No, estimator is given for \"" + legendre_fet_tally_name
            + "\", therefore, collision estimator will be taken, by default.");
    }
    std::string estimator_name = node["estimator"].as<std::string>();

    // Get the enrgy bounds, if any is given
    std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;
    if (node["energy-bounds"]){
        std::vector<double> energy_bounds = node["energy-bounds"].as<std::vector<double>>();
        energy_filter_ = std::make_shared<EnergyFilter>(energy_bounds);
    }

    // Get the cartesian type position filter
    std::shared_ptr<CartesianFilter> cartesian_filter_ 
                                    = make_cartesian_filter(node);


    // Get the Legendre-FET order
    if ( !node["FET-order"] ){
        fatal_error("Legendre-FET order is not given for the " 
            + legendre_fet_tally_name + "tally.");
    }
    int FET_order_int_type = node["FET-order"].as<int>();
    if (FET_order_int_type < 0 ){
        fatal_error("Legendre-FET order should not be negative.");
    }
    std::size_t FET_order_ = static_cast<std::size_t>(FET_order_int_type);

    // Get the Legendre-FET axes
    if ( !node["FET-axis"] ){
        fatal_error("The axes for " + legendre_fet_tally_name + " should be given. "); 
    }
    std::vector<std::string> fet_axis_type_string 
                            = node["FET-axis"].as<std::vector<std::string>>();
    
    std::vector<LegendreFET::Axis> fet_axis;
    fet_axis.reserve( fet_axis_type_string.size() );

    for ( auto& c : fet_axis_type_string ){
        if ( c == "X" )
            fet_axis.push_back(LegendreFET::Axis::X);

        if ( c == "Y" )
            fet_axis.push_back(LegendreFET::Axis::Y);

        if ( c == "Z" )
            fet_axis.push_back(LegendreFET::Axis::Z);
    }

    if ( !(fet_axis.size() == fet_axis_type_string.size()) ){
        fatal_error("A non-allowable name of the axis is given for " + legendre_fet_tally_name + ".");
    }

    // Make the Legendre FET tally class
    std::shared_ptr<ITally> itally_legendre_fet_tally_ = nullptr;
    if ( estimator_name == "collision" ){
        itally_legendre_fet_tally_ = std::make_shared<LegendreFET>(FET_order_, LegendreFET::Quantity::Flux, 
                                        LegendreFET::Estimator::Collision, legendre_fet_tally_name,
                                        fet_axis, cartesian_filter_, energy_filter_);
    }

    if ( estimator_name == "track-length" ){
        /*itally_legendre_fet_tally_ = std::make_shared<LegendreFET>(FET_order_, LegendreFET::Quantity::Flux, 
                                        LegendreFET::Estimator::TrackLength, legendre_fet_tally_name,
                                        fet_axis, cartesian_filter_, energy_filter_);
        */
        fatal_error("The track-length for the legendre-FET is not avilable as asked in " + legendre_fet_tally_name + " tally. ");
    }

    if ( itally_legendre_fet_tally_ == nullptr )
        fatal_error("Incorrect \"estimator\" is given.");


    return itally_legendre_fet_tally_;
}