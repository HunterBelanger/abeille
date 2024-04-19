#include <tallies/legendre_fet.hpp>
#include <utils/mpi.hpp>

#include <tallies/itally.hpp>



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
    

    }

void LegendreFET::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat ){
    std::vector<size_t> index_position;
    std::size_t index_E;

    if ( energy_in_->get_index(p.E(), index_E) ){

        index_position = cartesian_filter_->get_indices(tktr) ;

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

            
        indexes_.insert(indexes_.begin(), index_E);
    

        const size_t axis_index = indexes_.size();
        const size_t FET_index = indexes_.size() + 1;
        
        indexes_.push_back(axis_index);
        indexes_.push_back(FET_index);
        
        double beta_n;

        size_t it_axis = 0;
        for ( auto& c : axes){  // Loop over the different axis and indexing is done using the it_axis
            for (size_t i = 0; i<FET_order+1; i++){    // loop over differnt FET order

                switch (c){

                case LegendreFET::Axis::X: {
                    const double xmin_ = cartesian_filter_->x_min();
                    const double xmax_ = cartesian_filter_->x_max();
                    const double x = 2 * (tktr.r().x() - xmin_) / (xmax_ - xmin_);

                    beta_n = collision_score * legendre(i ,x); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
                    
            #ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;
                    it_axis++;
                    }
                    break;
                
                case LegendreFET::Axis::Y: {
                    const double ymin_ = cartesian_filter_->y_min();
                    const double ymax_ = cartesian_filter_->y_max();
                    const double y = 2 * (tktr.r().y() - ymin_) / (ymax_ - ymin_);

                    beta_n = collision_score * legendre(i ,y); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
            #ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;
                    it_axis++;
                    }
                    break;
                
                case LegendreFET::Axis::Z: {
                    const double zmin_ = cartesian_filter_->z_min();
                    const double zmax_ = cartesian_filter_->z_max();
                    const double z = 2 * (tktr.r().z() - zmin_) / (zmax_ - zmin_);

                    beta_n = collision_score * legendre(i , z); // score for i-th order's basis function
                    indexes_[axis_index] = it_axis;
                    indexes_[FET_index] = i;
            #ifdef ABEILLE_USE_OMP
            #pragma omp atomic
            #endif
                    tally_gen_score(indexes_) += beta_n;
                    it_axis++;
                    }
                    break;
                
                }
            }
        }
    }
}

