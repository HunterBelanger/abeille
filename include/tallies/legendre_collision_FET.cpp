#include "legendre_collision_FET.hpp"
#include "utils/mpi.hpp"

#include "ITally.h"



CollisionEstimatorFET::CollisionEstimatorFET(size_t FET_order_, Quantity quantity_,
    std::vector<Axis> axes_,
    std::shared_ptr<PositionFilter> position_filter_, 
    std::shared_ptr<EnergyFilter> energy_in_,
    std::shared_ptr<EnergyFilter> energy_out_)
    :
    FET_order(FET_order_), ITally(FET_order_, quantity_, axes_, position_filter_, energy_in_, energy_out_)
    {
    
    //for legendre, check the size of axis vector should be between than 3.
    if ( axes.size() >= 1 && axes.size()<=3 )
        std::cout<<"Fatal error, the no. of given axis is not between 1 to 3.\n";
    else{
        
    }

    }

void CollisionEstimatorFET::score_particle(const Particle& p, const Tracker& tktr, MaterialHelper& mat ){
    std::array<int, 3> index_position;
    std::size_t index_E;

    if ( position_filter->get_index(tktr.r(), index_position) 
    && energy_in->get_index(p.E(), index_E) ){

        const double Et = mat.Et(p.E());

        double collision_score = 1.0 / Et;

        switch(quantity){
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
        
        // adding the scores
        
        double beta_n;
        size_t it_axis = 0;
        for ( auto& c : axes){  // Loop over the different axis
            for (size_t i = 0; i<FET_order+1; i++){    // loop over differnt FET order

                switch (c)
                {
                case (Axis::X):
                    const double xmin_ = position_filter->x_min();
                    const double xmax_ = position_filter->x_max();
                    const double x = 2 * (tktr.r().x() - xmin_) / (xmax_ - xmin_);

                    const double beta_n = collision_score * legendre(i ,x); // score for i-th order's basis function

                    score_gen(index_E, 
                            index_position[0], index_position[1], index_position[2],
                            it_axis, i ) += beta_n;
                    it_axis++;
                    break;
                
                case (Axis::Y):
                    const double ymin_ = position_filter->y_min();
                    const double ymax_ = position_filter->y_max();
                    const double y = 2 * (tktr.r().y() - ymin_) / (ymax_ - ymin_);

                    const double beta_n = collision_score * legendre(i ,y); // score for i-th order's basis function

                    score_gen(index_E, 
                            index_position[0], index_position[1], index_position[2],
                            it_axis, i ) += beta_n;
                    it_axis++;
                    break;
                
                case (Axis::Z):
                    const double zmin_ = position_filter->z_min();
                    const double zmax_ = position_filter->z_max();
                    const double z = 2 * (tktr.r().z() - zmin_) / (zmax_ - zmin_);

                    const double beta_n = collision_score * legendre(i ,z); // score for i-th order's basis function

                    score_gen(index_E, 
                            index_position[0], index_position[1], index_position[2],
                            it_axis, i ) += beta_n;
                    it_axis++;
                    break;
                
                default:
                    break;
                }
            }
        }
    }
}


void CollisionEstimatorFET::record_variance(double multiplier){
        gen++;
        const double dg = static_cast<double>(gen);
        const double invs_dg = 1. / dg;

        // All worker threads must send their generation score to the master.
        // Master must recieve all generations scores from workers and add
        // them to it's own generation score.
        mpi::Reduce_sum(score_gen.data_vector(), 0);

        // Only try to update average and variance is we are master, as worker
        // processes don't have copies of this data, so it will seg-fault.
        if (mpi::rank == 0) {
        #ifdef ABEILLE_USE_OMP
        #pragma omp parallel for schedule(static)
        #endif
            for (size_t i = 0; i < FET_coeff.size(); i++) {
                // Get new average
                double old_avg = FET_coeff[i];
                double val = score_gen[i] * multiplier;
                double avg = old_avg + (val - old_avg) * invs_dg;
                FET_coeff[i] = avg;

                // Get new variance
                double var = FET_var[i];
                var = var + ((val - old_avg) * (val - avg) - (var)) * invs_dg;
                FET_var[i] = var;
                }
        }  
    }