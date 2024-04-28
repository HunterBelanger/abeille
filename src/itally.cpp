#include <tallies/itally.hpp>
#include <utils/output.hpp>

#include <tallies/tallies.hpp>
#include <tallies/general_tally.hpp>
#include <tallies/legendre_fet.hpp>
#include <tallies/box_position_filter.hpp>
#include <tallies/mesh_position_filter.hpp>

//#include <boost/container/static_vector.hpp>

void ITally::record_generation(double multiplier){
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
            //std::cout<<"In record tally "<<tally_avg[i]<<"\n";
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
        return "track-length";
        break;

    default:
        return "unkown";
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

void make_itally(Tallies& tallies, const YAML::Node& node){

    if ( !node["name"] ){
        fatal_error("Tally name is not given.");
    }
    std::string tally_name_ = node["name"].as<std::string>();
    std::string tally_type_ = "general";
    if (!node["tally-type"]){
        warning("The tally-type is not given therefore, \"general\" will be assumed.");
    }
    tally_type_ = node["tally-type"].as<std::string>();

    std::shared_ptr<ITally> new_ITally = nullptr;
    if ( tally_type_ == "general" ){
        new_ITally =  make_general_tally(node);
    }

    if ( tally_type_ == "legendre-FET" ){
        new_ITally = make_legendre_fet(node);
    }

    if ( tally_type_ == "zernike-FET" ){
        fatal_error("The zernike-FET for " + tally_name_ + " is not supported yet.");
    }

    if ( new_ITally == nullptr ){
        fatal_error("For the " + tally_name_ + ", something went worng in the input.");
    }

    // Add the new_ITally of type ITally into the "tallies"
    tallies.add_ITally(new_ITally);

}

/*
void make_itally(Tallies& tallies, const YAML::Node& node){
    std::string tally_name_;
    if( !node["name"] ){
        fatal_error("No valid name is provided.");
    }
    tally_name_ = node["name"].as<std::string>();

    std::string tally_type = "general";

    size_t tally_fet_order = 0;

    std::vector<std::string> legendre_fet_vec;
    legendre_fet_vec.reserve(3);

    std::vector<LegendreFET::Axis> axes_vec;
    axes_vec.reserve(3);
    
    if (node["tally-type"]){
        tally_type = node["tally-type"].as<std::string>();
        // Check for the correct input
        if ( !(tally_type == "general" || 
            tally_type == "legendre-FET" || tally_type == "zernike-FET") ){
            fatal_error("No valid tally-type is provided.");
        }
    }

    if ( tally_type == "legendre-FET" || tally_type == "zernike-FET"){
        if (node["FET-order"]){
              tally_fet_order = node["FET-order"].as<std::size_t>();
        }
        else{    fatal_error("FET order is not provided.");
        }

        if (node["FET-axis"]){
            legendre_fet_vec = node["FET-axis"].as<std::vector<std::string>>();
            bool check_axes = false;
            for(auto& c_ : legendre_fet_vec){
                if ( c_ == "X"){
                    axes_vec.push_back(LegendreFET::Axis::X);
                    check_axes = true;
                }
                
                if ( c_ == "Y"){
                    axes_vec.push_back(LegendreFET::Axis::Y);
                    check_axes = true;
                }
                
                if ( c_ == "Z"){
                    axes_vec.push_back(LegendreFET::Axis::Z);
                    check_axes = true;
                }
                    
            }
            
            if (check_axes == false)    fatal_error("FET evaluation axis is not provided.");

        }
        else    fatal_error("FET evaluation axis is not provided.");

    }

    std::string quant;
    
    if (!node["quantity"])
        fatal_error("No, quantity is given for evaluation.");
    
    quant = node["quantity"].as<std::string>();


    std::string estimator_ = "collision";
    if (!node["estimator"])
        warning("No, estimator is given for \""+tally_name_
                +"\", therefore, collision estimator will be taken, by default.");
    else
        estimator_ = node["estimator"].as<std::string>();
    
    std::string pos_filter_name= "box";
    Position r_low_, r_high_;
    std::vector<size_t> pos_dime;
    
    if( node["Position-Filter"]){
        pos_filter_name =  node["Position-Filter"].as<std::string>();

        if (pos_filter_name == "box" || pos_filter_name == "regular-reactangular-mesh"){
        std::vector<double> loc = node["low"].as<std::vector<double>>();
        // if (loc.size() != 3) { fatal_error("Low has not three elements.")}
        r_low_ = Position(loc[0], loc[1], loc[2]);

        std::vector<double> hi = node["high"].as<std::vector<double>>();
        // if (hi.size() != 3) { fatal_error("High has not three elements.")}
        r_high_ = Position( hi[0], hi[1], hi[2]);
        }
        if(node["shape"]){
            pos_dime = node["shape"].as<std::vector<size_t>>();
        }  
    }

    std::shared_ptr<PositionFilter> pos_filter;
    std::shared_ptr<CartesianFilter> cart_filter;

    /*if (pos_filter_name == "box" ){
        pos_filter = std::make_shared<BoxPositionFilter>(r_low_, r_high_);
        cart_filter = std::make_shared<BoxPositionFilter>(r_low_, r_high_);
    }else if (pos_filter_name == "reactangular-mesh"){
        pos_filter = std::make_shared<MeshPositionFilter>(r_low_, r_high_, pos_dime[0], pos_dime[1], pos_dime[2]);
        cart_filter = std::make_shared<MeshPositionFilter>(r_low_, r_high_, pos_dime[0], pos_dime[1], pos_dime[2]);
    }

    pos_filter = make_position_filter(node);


    std::vector<double> ebounds;
    if (node["energy-bounds"]){
        ebounds = node["energy-bounds"].as<std::vector<double>>();
    }

    std::shared_ptr<EnergyFilter> energy_filter_ = std::make_shared<EnergyFilter>(ebounds);

    std::shared_ptr<ITally> new_tally;
    if ( estimator_ == "collision" ){
       
        if(tally_type == "general"){
                //new_tally = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux, GeneralTally::Estimator::Collision, tally_name_,
                //                        pos_filter, energy_filter_);
                new_tally = make_general_tally(node);
                std::cout<<"weis made.\n";
        }

        if(tally_type == "legendre-FET"){
                cart_filter = make_cartesian_filter(node);
                /*new_tally = std::make_shared<LegendreFET>(tally_fet_order, LegendreFET::Quantity::Flux, LegendreFET::Estimator::Collision, tally_name_,
                                        axes_vec, cart_filter, energy_filter_);
                
                new_tally = make_legendre_fet(node);
        }

    }

    if ( estimator_ == "track-length" ){
        new_tally 
            = std::make_shared<GeneralTally>(GeneralTally::Quantity::Flux, GeneralTally::Estimator::TrackLength, tally_name_,
                                             pos_filter, energy_filter_);
    }
      
    tallies.add_ITally(new_tally);
    std::cout<<">>>>>>> "<< new_tally->estimator_str()<<"\n";
    
    
}
*/