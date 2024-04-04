#include <tallies/itally.hpp>
#include <utils/output.hpp>

//#include <boost/container/static_vector.hpp>

void ITally::record_generation(double multiplier){
    gen_++;
    std::cout<<"Itally gen = "<<gen_<<"\n";
    
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
            std::cout<<"In record tally "<<tally_avg[i]<<"\n";
            }
    }
        
    // Clear the entry for the tally_gen
    tally_gen_score.fill(0.0);    
}

std::string ITally::estimator_str(){
    switch (estimator_)
    {
    case Estimator::Collision:
        return "collision";
        break;
    
    case Estimator::TrackLength:
        return "tracklength";
        break;

    default:
        break;
    }
}

void ITally::write_tally(){
    // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + this->tally_name);

  // First write coordinates and number of groups
  /*std::vector<double> x_bounds(Nx + 1, 0.);
  for (std::size_t i = 0; i <= Nx; i++) {
    x_bounds[i] = (static_cast<double>(i) * dx) + r_low.x();
  }
  tally_grp.createAttribute("x-bounds", x_bounds);

  std::vector<double> y_bounds(Ny + 1, 0.);
  for (std::size_t i = 0; i <= Ny; i++) {
    y_bounds[i] = (static_cast<double>(i) * dy) + r_low.y();
  }
  tally_grp.createAttribute("y-bounds", y_bounds);

  std::vector<double> z_bounds(Nz + 1, 0.);
  for (std::size_t i = 0; i <= Nz; i++) {
    z_bounds[i] = (static_cast<double>(i) * dz) + r_low.z();
  }
  tally_grp.createAttribute("z-bounds", z_bounds);
    
   
  tally_grp.createAttribute("energy-bounds", energy_bounds);
 
  // Save the quantity
  tally_grp.createAttribute("quantity", this->quantity_str());

  if (this->quantity_str() == "mt") {
    tally_grp.createAttribute("mt", this->mt());
  }*/

  // Save the estimator
  tally_grp.createAttribute("estimator", this->estimator_str());

  // Convert flux_var to the error on the mean
  for (size_t l = 0; l < tally_var.size(); l++)
    tally_var[l] = std::sqrt(tally_var[l] / static_cast<double>(gen_));

  // Add data sets for the average and the standard deviation
  auto avg_dset =
      tally_grp.createDataSet<double>("avg", H5::DataSpace(tally_avg.shape()));
  avg_dset.write_raw(&tally_avg[0]);

  auto std_dset =
      tally_grp.createDataSet<double>("std", H5::DataSpace(tally_var.shape()));
  std_dset.write_raw(&tally_var[0]);

}