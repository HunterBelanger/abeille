#include <tallies/ITally.h>

void ITally::record_tally(double multiplier){
    gen_++;
    const double dg = static_cast<double>(gen_);
    const double invs_dg = 1. / dg;

    // All worker threads must send their generation score to the master.
    // Master must recieve all generations scores from workers and add
    // them to it's own generation score.
    mpi::Reduce_sum(tally_gen_score.data_vector(), 0);

    // Only try to update average and variance is we are master, as worker
    // processes don't have copies of this data, so it will seg-fault.
    if (mpi::rank == 0) {
    #ifdef ABEILLE_USE_OMP
    #pragma omp parallel for schedule(static)
    #endif
        for (size_t i = 0; i < tally_avg.size(); i++) {
            // Get new average

            double old_avg = tally_avg[i];
            double val = tally_gen_score[i] * multiplier;
            double avg = old_avg + (val - old_avg) * invs_dg;
            tally_avg[i] = avg;

            // Get new variance
            double var = tally_var[i];
            var = var + ((val - old_avg) * (val - avg) - (var)) * invs_dg;
            tally_var[i] = var;
            std::cout<<"In record tally "<<tally_avg[i]<<","<<tally_avg.size()<<"\n";
            }
    }
        
    // Clear the entry for the tally_gen
    tally_gen_score.fill(0.0);    
}
