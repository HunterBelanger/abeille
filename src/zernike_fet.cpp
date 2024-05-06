#include <tallies/zernike_fet.hpp>
#include <utils/error.hpp>
#include <utils/legendre.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<std::size_t, 6>;

ZernikeFET::ZernikeFET(std::shared_ptr<CylinderPositionFilter> cylinder_filter,
                    std::shared_ptr<EnergyFilter> energy_filter,
                    std::size_t zernike_order, std::size_t lengendre_order,
                    ZernikeFET::Quantity quantity, ZernikeFET::Estimator estimator, std::string name)
    :
    ITally(quantity, estimator, name),
    cylinder_filter_(cylinder_filter),
    energy_filter_(energy_filter),
    zr_polynomial_(),
    zr_order_(zernike_order),
    legen_order_(lengendre_order),
    axial_direction_(){
    StaticVector6 tally_shape;
    // add the dimension for energy_in_ only if exist
    if (energy_filter_) {
        std::size_t ne = energy_filter_->size();
        tally_shape.push_back(ne);
    }

    // get the shape or dimensions for the cylinder-filter
    // the zernike-fet requires the co-ordinate information
    if (cylinder_filter_ == nullptr) {
    fatal_error("for tally " + name + ", the cylinder-filte's information is required for zernike-fet.");
    }
    StaticVector3 cylinder_shape = cylinder_filter_->get_shape();
    tally_shape.insert(tally_shape.end(), cylinder_shape.begin(), cylinder_shape.end());

    // another dimnesion of length two, which is related to the zernike responsible 
    // for radial and legendre reposible for the axis. first will be zernike and 
    // second will be legendre.
    tally_shape.push_back(2);

    // last-dimension is for the orders, however, the radial will incorporated by
    // zernike and axial will be legendre, therefore, for both the order can be 
    // same or different, therefore, a max of both will be taken, so all the tally
    // will remain at one place.  
    std::size_t max_order = zr_order_;
    if ( legen_order_ > zr_order_ ){
        max_order = legen_order_;
    }
    tally_shape.push_back(max_order+1);

    // reallocate and fill with zeros for the taly avg, gen-score and varaince
    tally_avg_.reallocate(tally_shape);
    tally_avg_.fill(0.0);

    tally_gen_score_.reallocate(tally_shape);
    tally_gen_score_.fill(0.0);

    tally_var_.reallocate(tally_shape);
    tally_var_.fill(0.0);

    // construct the zernike polynomials upto the given order
    zr_polynomial_ = ZernikePolynomials(zr_order_);
}

void ZernikeFET::score_collision(const Particle& p, const Tracker& tktr, MaterialHelper& mat){
    StaticVector6 indices;
    // get the energy-index, if energy-filter exists
    if (energy_filter_) {
        std::optional<std::size_t> E_indx = energy_filter_->get_index(p.E());

        if (E_indx.has_value()) {
        std::size_t index_E = E_indx.value();
        indices.push_back(index_E);
        } else
        // Not inside any energy bin. Don't score.
        return;
    }

    // get the positional indexes from cylinder_filter 
    StaticVector3 cylinder_index = cylinder_filter_->get_indices(tktr);
    if (cylinder_index.empty()){
        // Not inside any energy bin. Don't score.
        return;
    }
    indices.insert(indices.begin(), cylinder_index.begin(), cylinder_index.end());

    // add one dimension for the different polynomial
    std::size_t poly_index = indices.size();
    indices.push_back(0);

    // add one dimension for the different order
    std::size_t FET_index = indices.size();
    indices.push_back(0);

    // get the particle base score for the collision
    const double Et = mat.Et(p.E());
    const double collision_score = particle_base_score(p, mat) / Et;

    // varaible for scoring
    double beta_n;

    // first- score for the zernike 
    std::pair<double, double> scaled_r_and_theta = cylinder_filter_->get_scaled_radius_and_angle(cylinder_index, tktr.r());
    const double scaled_r = scaled_r_and_theta.first;
    const double theta = scaled_r_and_theta.second;
    const std::vector<double> zr_value = zr_polynomial_.evaluate_zernikes(scaled_r, theta);
    if ( zr_value.size() != (zr_order_+1)){
        fatal_error("for the zernike-FET, something happened as the sizes for evaluated polynomials and zr-order+1 is not same.");
    }
    for(std::size_t i = 0; i<=zr_order_; i++){
        // score for i-th order's basis function
        beta_n = collision_score * zr_value[i];
        indices[FET_index] = i;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_(indices) += beta_n;
    }

    // second- score for the legendre
    indices[poly_index] = 1;
    const double zmin = cylinder_filter_->z_min(cylinder_index);
    const double zmax = cylinder_filter_->z_max(cylinder_index);
    double z = tktr.r().z();
    if ( axial_direction_ == CylinderPositionFilter::Orientation::Y){
        z = tktr.r().y();
    } else if ( axial_direction_ == CylinderPositionFilter::Orientation::X){
        z = tktr.r().z();
    }
    const double scaled_z = 2 * ( z - zmin) / ( zmin-zmax) - 1.0; 
    for( std::size_t i = 0; i<=legen_order_; i++){
        // score for i-th order's basis function
        beta_n = collision_score * legendre(i, scaled_z);
        indices[FET_index] = i;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_(indices) += beta_n;    
    }

}

void ZernikeFET::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + tally_name_);
  /*
    // First write coordinates and number of groups
    std::vector<double> x_bounds(Nx + 1, 0.);
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
  */
  // Save the quantity
  tally_grp.createAttribute("quantity", quantity_str());

  /*if (this->quantity_str() == "mt") {
    tally_grp.createAttribute("mt", this->mt());
  }*/

  // Save the estimator
  tally_grp.createAttribute("estimator", estimator_str());

  // Convert flux_var to the error on the mean
  for (size_t l = 0; l < tally_var_.size(); l++)
    tally_var_[l] = std::sqrt(tally_var_[l] / static_cast<double>(gen_));

  // Add data sets for the average and the standard deviation
  auto avg_dset =
      tally_grp.createDataSet<double>("avg", H5::DataSpace(tally_avg_.shape()));
  avg_dset.write_raw(&tally_avg_[0]);

  auto std_dset =
      tally_grp.createDataSet<double>("std", H5::DataSpace(tally_var_.shape()));
  std_dset.write_raw(&tally_var_[0]);
}

// make the tally-zernike
std::shared_ptr<ZernikeFET> make_zernike_fet(const YAML::Node& node){
    // Check the name of the tallt is given or not.
    if (!node["name"] && !node["name"].IsScalar()) {
        fatal_error("No valid name is provided.");
    }
    std::string zernike_fet_tally_name = node["name"].as<std::string>();

    // Check the estimator is given or not.
    if (!node["estimator"] && !node["estimator"].IsScalar()) {
        fatal_error("No valid estimator is given for " + zernike_fet_tally_name +
                    " tally.");
    }
    std::string estimator_name = node["estimator"].as<std::string>();

    // Check for the quantity
    std::string given_quantity = "";
    if (!node["quantity"] && !node["quantity"].IsScalar()) {
        fatal_error("No valid quantity is given for " + zernike_fet_tally_name +
                    " tally.");
    }
    given_quantity = node["quantity"].as<std::string>();
    ZernikeFET::Quantity quant;

    if (given_quantity == "flux") {
        quant = ZernikeFET::Quantity::Flux;
    } else if (given_quantity == "fission") {
        quant = ZernikeFET::Quantity::Fission;
    } else if (given_quantity == "absorption") {
        quant = ZernikeFET::Quantity::Absorption;
    } else if (given_quantity == "elastic") {
        quant = ZernikeFET::Quantity::Elastic;
    } else {
        fatal_error("For " + zernike_fet_tally_name +
                    " tally, a unknown quantity is given.");
    }

    // Get the enrgy filter, if any is given
    std::shared_ptr<EnergyFilter> energy_filter = nullptr;
    if (node["energy-bounds"]) {
        if (!node["energy-bounds"].IsSequence()) {
        fatal_error("energy-bounds is not provided in the list format.");
        }
        energy_filter = make_energy_filter(node);
    }

    // Get the cylinder type position filter
    if (!node["position-filter-type"] || !node["position-filter-type"].IsScalar()) {
        fatal_error("For " + zernike_fet_tally_name +
                    ", a valid position-filter must be given.");
    }
    std::shared_ptr<CylinderPositionFilter> cylinder_filter =
        make_cylinder_position_filter(node);

    // get the order of zernike 
    if (!node["zernike-order"] && !node["zernike-order"].IsScalar()){
        fatal_error("a valid zernike-order is not given for " +
                zernike_fet_tally_name + "tally.");
    }
    std::size_t zernike_fet_order = node["zernike-order"].as<std::size_t>();

    // get the order of legendre
    if (!node["legendre-order"] && !node["legendre-order"].IsScalar()){
        fatal_error("a valid legendre-order is not given for " +
                zernike_fet_tally_name + "tally.");
    }
    std::size_t legendre_fet_order = node["legendre-order"].as<std::size_t>();

    // make the ZernikeFET class for the given tally description
    std::shared_ptr<ZernikeFET> itally_zernike_fet_tally = nullptr;
    if (estimator_name == "collision" ){
        itally_zernike_fet_tally = std::make_shared<ZernikeFET>(cylinder_filter, energy_filter, 
                                zernike_fet_order, legendre_fet_order, quant, 
                                ZernikeFET::Estimator::Collision,  zernike_fet_tally_name);
    } else if ( estimator_name == "track-length"){
            fatal_error(
        "The track-length for the zernike-fet is not supported as asked in " +
        zernike_fet_tally_name + " tally. ");
    } else {
        fatal_error("Incorrect \"estimator\" is given for the " +
                zernike_fet_tally_name + " tally.");
    }

    return itally_zernike_fet_tally;
}