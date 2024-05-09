#include <tallies/tallies.hpp>
#include <tallies/zernike_fet.hpp>
#include <utils/error.hpp>
#include <utils/legendre.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<std::size_t, 6>;

ZernikeFET::ZernikeFET(std::shared_ptr<CylinderFilter> cylinder_filter,
                       std::shared_ptr<EnergyFilter> energy_filter,
                       std::size_t zernike_order, std::size_t lengendre_order,
                       Quantity quantity, Estimator estimator, std::string name)
    : ITally(quantity, estimator, name),
      cylinder_filter_(cylinder_filter),
      energy_filter_(energy_filter),
      zr_polynomial_(zernike_order),
      zr_order_(zernike_order),
      legen_order_(lengendre_order),
      axial_direction_() {
  StaticVector6 tally_shape;
  // add the dimension for energy_in_ only if exist
  if (energy_filter_) {
    std::size_t ne = energy_filter_->size();
    tally_shape.push_back(ne);
  }

  // get the shape or dimensions for the cylinder-filter
  // the zernike-fet requires the co-ordinate information
  if (cylinder_filter_ == nullptr) {
    fatal_error("ZernikeFET has nullptr cylinder-filer.");
  }
  StaticVector3 cylinder_shape = cylinder_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), cylinder_shape.begin(),
                     cylinder_shape.end());
  axial_direction_ = cylinder_filter_->get_axial_direction();

  // another dimnesion of length two, which is related to the zernike
  // responsible for radial and legendre reposible for the axis. first will be
  // zernike and second will be legendre.
  tally_shape.push_back(2);

  // last-dimension is for the orders, however, the radial will incorporated by
  // zernike and axial will be legendre, therefore, for both the order can be
  // same or different, therefore, a max of both will be taken, so all the tally
  // will remain at one place.
  std::size_t max_order = zr_order_;
  if (legen_order_ > zr_order_) {
    max_order = legen_order_;
  }
  tally_shape.push_back(max_order + 1);

  // reallocate and fill with zeros for the taly avg, gen-score and varaince
  tally_avg_.reallocate(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.reallocate(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.reallocate(tally_shape);
  tally_var_.fill(0.0);
}

ZernikeFET::ZernikeFET(std::shared_ptr<CylinderFilter> cylinder_filter,
                       std::shared_ptr<EnergyFilter> energy_filter,
                       std::size_t zernike_order, Quantity quantity,
                       Estimator estimator, std::string name)
    : ITally(quantity, estimator, name),
      cylinder_filter_(cylinder_filter),
      energy_filter_(energy_filter),
      zr_polynomial_(zernike_order),
      zr_order_(zernike_order),
      legen_order_(),
      axial_direction_() {
  StaticVector6 tally_shape;
  // add the dimension for energy_in_ only if exist
  if (energy_filter_) {
    std::size_t ne = energy_filter_->size();
    tally_shape.push_back(ne);
  }

  // get the shape or dimensions for the cylinder-filter
  // the zernike-fet requires the co-ordinate information
  if (cylinder_filter_ == nullptr) {
    fatal_error("ZernikeFET has nullptr cylinder-filer.");
  }
  StaticVector3 cylinder_shape = cylinder_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), cylinder_shape.begin(),
                     cylinder_shape.end());
  axial_direction_ = cylinder_filter_->get_axial_direction();
  // another dimnesion of length two, which is related to the zernike
  // responsible for radial and legendre reposible for the axis. first will be
  // zernike and second will be legendre.
  tally_shape.push_back(1);

  // last-dimension is for the orders, however, the radial will incorporated by
  // zernike and axial will be legendre, therefore, for both the order can be
  // same or different, therefore, a max of both will be taken, so all the tally
  // will remain at one place.
  std::size_t max_order = zr_order_;
  tally_shape.push_back(max_order + 1);

  // in-this particular constructor, there will no legendre-order
  check_for_legendre = false;

  // reallocate and fill with zeros for the taly avg, gen-score and varaince
  tally_avg_.reallocate(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.reallocate(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.reallocate(tally_shape);
  tally_var_.fill(0.0);
}

void ZernikeFET::score_collision(const Particle& p, const Tracker& tktr,
                                 MaterialHelper& mat) {
  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_filter_) {
    std::optional<std::size_t> E_indx = energy_filter_->get_index(p.E());
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    indices.push_back(E_indx.value());
  }

  // get the positional indexes from cylinder_filter
  StaticVector3 cylinder_index = cylinder_filter_->get_indices(tktr);
  if (cylinder_index.empty()) {
    // Not inside any energy bin. Don't score.
    return;
  }
  indices.insert(indices.end(), cylinder_index.begin(), cylinder_index.end());

  // add one dimension for the different polynomial
  std::size_t poly_index = indices.size();
  indices.push_back(0);

  // add one dimension for the different order
  std::size_t FET_index = indices.size();
  indices.push_back(0);

  // Get the particle base score for the collision
  const double Et = mat.Et(p.E());
  const double collision_score = particle_base_score(p, mat) / Et;

  // varaible for scoring
  double beta_n;

  // first- score for the zernike
  std::pair<double, double> scaled_r_and_theta =
      cylinder_filter_->get_scaled_radius_and_angle(cylinder_index, tktr.r());
  const double scaled_r = scaled_r_and_theta.first;
  const double theta = scaled_r_and_theta.second;
  const std::vector<double> zr_value =
      zr_polynomial_.evaluate_zernikes(scaled_r, theta);

  for (std::size_t i = 0; i <= zr_order_; i++) {
    // score for i-th order's basis function
    beta_n = collision_score * zr_value[i];
    indices[FET_index] = i;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_(indices) += beta_n;
  }

  // second- score for the legendre
  if (check_for_legendre == true) {
    indices[poly_index] = 1;
    const double zmin = cylinder_filter_->z_min(cylinder_index);
    const double zmax = cylinder_filter_->z_max(cylinder_index);
    double z = tktr.r().z();
    if (axial_direction_ == CylinderFilter::Orientation::Y) {
      z = tktr.r().y();
    } else if (axial_direction_ == CylinderFilter::Orientation::X) {
      z = tktr.r().z();
    }
    const double scaled_z = 2 * (z - zmin) / (zmin - zmax) - 1.0;
    for (std::size_t i = 0; i <= legen_order_; i++) {
      // score for i-th order's basis function
      beta_n = collision_score * legendre(i, scaled_z);
      indices[FET_index] = i;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_(indices) += beta_n;
    }
  }
}

void ZernikeFET::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + tally_name_);

  // Save the type
  tally_grp.createAttribute("type", "zernike-fet");

  // Save the quantity
  tally_grp.createAttribute("quantity", quantity_str());

  // Save the estimator
  tally_grp.createAttribute("estimator", estimator_str());

  // Save the energy-in filter
  if (energy_filter_) {
    tally_grp.createAttribute("energy-filter", energy_filter_->id());
  }

  // Save position filter
  tally_grp.createAttribute("position-filter", cylinder_filter_->id());

  // Save FET order
  tally_grp.createAttribute("zernike-order", zr_order_);

  if (check_for_legendre)
    tally_grp.createAttribute("legendre-order", legen_order_);

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
std::shared_ptr<ZernikeFET> make_zernike_fet(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || node["name"].IsScalar() == false) {
    fatal_error("No valid name is provided on tally.");
  }
  std::string name = node["name"].as<std::string>();

  // Check the estimator is given or not. We default to collision estimators.
  std::string estimator_name = "collision";
  if (node["estimator"] && node["estimator"].IsScalar()) {
    estimator_name = node["estimator"].as<std::string>();
  } else if (node["estimator"]) {
    fatal_error("Invalid estimator entry is given on tally \"" + name + "\".");
  }
  Estimator estimator;
  if (estimator_name == "collision") {
    estimator = Estimator::Collision;
  } else if (estimator_name == "track-length") {
    estimator = Estimator::TrackLength;
    fatal_error(
        "On tally " + name +
        ", track-length estimator is not yet supported on zernike-fet.");
  } else {
    std::stringstream mssg;
    mssg << "The tally " << name << " was provided with an unkown estimator \""
         << estimator_name << "\".";
    fatal_error(mssg.str());
  }

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    fatal_error("No quantity is given for tally " + name + ".");
  }
  given_quantity = node["quantity"].as<std::string>();
  Quantity quant = read_quantity(given_quantity, name);

  // Get the tallies instance
  auto& tallies = Tallies::instance();

  // Get the energy filter, if any is given
  std::shared_ptr<EnergyFilter> energy_filter = nullptr;
  if (node["energy-filter"] && node["energy-filter"].IsScalar()) {
    std::size_t energy_id = node["energy-filter"].as<std::size_t>();
    energy_filter = tallies.get_energy_filter(energy_id);
    if (energy_filter == nullptr) {
      std::stringstream mssg;
      mssg << "For tally " << name << ", cannot find energy filter with id "
           << energy_id << ".";
      fatal_error(mssg.str());
    }
  } else if (node["energy-filter"]) {
    fatal_error("Invalid energy-filter entry on tally " + name + ".");
  }

  // Get the cylinder type position filter
  std::shared_ptr<CylinderFilter> cylinder_filter = nullptr;
  if (!node["position-filter"] || node["position-filter"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Zernike-FET tally " << name
         << " has invalid/nonexistent position-filter entry.";
    fatal_error(mssg.str());
  }
  std::size_t position_id = node["position-filter"].as<std::size_t>();
  cylinder_filter = tallies.get_cylinder_filter(position_id);
  if (cylinder_filter == nullptr) {
    std::stringstream mssg;
    mssg << "Zernike-FET tally " << name << " was provided position-filter id "
         << position_id << ". No cylinder filter with this id was found.";
    fatal_error(mssg.str());
  }

  // Get the Zernike order
  if (!node["zernike-order"] || node["zernike-order"].IsScalar() == false) {
    fatal_error("a valid zernike-order is not given for " + name + "tally.");
  }
  std::size_t zernike_fet_order = node["zernike-order"].as<std::size_t>();

  // get the order of legendre, if exists
  bool legendre_exist = true;
  std::size_t legendre_fet_order;
  if (!node["legendre-order"]) {
    legendre_exist = false;
  } else if (!node["legendre-order"].IsScalar()) {
    std::stringstream mssg;
    mssg << "Tally " << name << " has invalid legendre-order entry.";
    fatal_error(mssg.str());
  } else {
    legendre_fet_order = node["legendre-order"].as<std::size_t>();
  }

  // If we have both legendre and zernike then use the first constructor,
  if (legendre_exist) {
    return std::make_shared<ZernikeFET>(cylinder_filter, energy_filter,
                                        zernike_fet_order, legendre_fet_order,
                                        quant, estimator, name);
  } else {
    // if we have the zernike only, then use the second constructor
    return std::make_shared<ZernikeFET>(cylinder_filter, energy_filter,
                                        zernike_fet_order, quant, estimator,
                                        name);
  }
}