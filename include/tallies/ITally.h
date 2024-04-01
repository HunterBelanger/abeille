#ifndef ITALLY_H
#define ITALLY_H

#include <memory>

#include <position_filter.h>
#include <energy_filter.h>
#include <ndarray.hpp>

class ITally{
    public: 
    enum class Quantity{
        Flux,
        Fission,
        Absorption,
        Elastic, 
    };

    enum class Axis{
        X, Y, Z
    };


    ITally(Quantity quantity_,
        std::shared_ptr<EnergyFilter> energy_in_,
        std::shared_ptr<EnergyFilter> energy_out_ = nullptr) 
        :
        FET_coeff(), FET_var(), score_gen(), quantity(quantity_),
        axes(axes_), position_filter(position_filter_), 
        energy_in(energy_in_), energy_out(energy_out_)
        {
            //for legendre, check the size of axis vector should be between than 3.
            if ( axes.size() >= 1 && axes.size()<=3 )
                std::cout<<"Fatal error, the no. of given axis is not between 1 to 3.\n";
            else{
                
            }
            
            uint64_t Nx = static_cast<uint64_t> (position_filter->get_Nx());
            uint64_t Ny = static_cast<uint64_t> (position_filter->get_Ny());
            uint64_t Nz = static_cast<uint64_t> (position_filter->get_Nz());
            size_t ND  = axes.size();
            uint64_t NE = static_cast<uint64_t> (energy_in->get_size());
            // Order+1 is required as these many element is required. 
            size_t N_FET_element = static_cast<size_t> (FET_order_ + 1);

            FET_coeff.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_coeff.fill(0.0);

            score_gen.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            score_gen.fill(0.0);

            FET_var.reallocate({NE, Nx, Ny, Nz, ND, N_FET_element});
            FET_var.fill(0.0);
            //FET_std.reallocate({NE, Nx, Ny, Nz});
            //FET_std.fill( UpperTriangularMatrix(FET_order + 1) );

            
        }
    
    ~ITally() = default;

    //For collision estmator
    virtual score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat );

    //For track-length estimator
    virtual score_flight(const Particle& p, const Tracker& tktr, MaterialHelper& mat );

    //Record the avg and variance for the generation
    void record_tally(double mulitplier);

     protected:

        std::shared_ptr<EnergyFilter> energy_in {nullptr};
        std::shared_ptr<EnergyFilter> energy_out {nullptr};
        
        NDArray<double> tally_avg;
        NDArray<double> talley_gen_score;
        NDArray<double> tally_var;

        size_t gen_ = 0;

        Quantity quantity;
        
        
};

void ITally::record_tally(double multiplier){
    gen_++;
    const double dg = static_cast<double>(gen_);
    const double invs_dg = 1. / dg;

    // All worker threads must send their generation score to the master.
    // Master must recieve all generations scores from workers and add
    // them to it's own generation score.
    mpi::Reduce_sum(talley_gen_score.data_vector(), 0);

    // Only try to update average and variance is we are master, as worker
    // processes don't have copies of this data, so it will seg-fault.
    if (mpi::rank == 0) {
    #ifdef ABEILLE_USE_OMP
    #pragma omp parallel for schedule(static)
    #endif
        for (size_t i = 0; i < tally_avg.size(); i++) {
            // Get new average
            double old_avg = tally_avg[i];
            double val = talley_gen_score[i] * multiplier;
            double avg = old_avg + (val - old_avg) * invs_dg;
            tally_avg[i] = avg;

            // Get new variance
            double var = tally_var[i];
            var = var + ((val - old_avg) * (val - avg) - (var)) * invs_dg;
            tally_var[i] = var;
            }
    }

    // Clear the entry for the tally_gen
    talley_gen_score.fill(0.0);    
}

#endif