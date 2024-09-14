#include <tallies/tallies.hpp>
#include <tallies/zernike_fet.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<std::size_t, 6>;

ZernikeFET::ZernikeFET(std::shared_ptr<CylinderFilter> cylinder_filter,
                       std::shared_ptr<EnergyFilter> energy_filter,
                       std::size_t zernike_order, std::size_t legendre_order,
                       Quantity quantity, Estimator estimator, std::string name)
    : ITally(quantity, estimator, name),
      cylinder_filter_(cylinder_filter),
      energy_filter_(energy_filter),
      zr_polynomial_(zernike_order),
      zr_order_(zernike_order),
      legen_order_(legendre_order),
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

  if ((quantity_.type == Quantity::Type::Source ||
       quantity_.type == Quantity::Type::RealSource ||
       quantity_.type == Quantity::Type::ImagSource) &&
      estimator_ != Estimator::Source) {
    std::stringstream mssg;
    mssg << "Tally " << tally_name_
         << " has a source-like qantity but does not use a source estimator.";
    fatal_error(mssg.str());
  }

  // throw the error if the infinte-cylinder and legendre-fet is given
  if (cylinder_filter_->is_infinite_cylinder() == true) {
    std::stringstream mssg;
    mssg << "Tally " << name
         << "cannot evaluate the legendre-fet over a infinte-cylinder with "
            "position-filter id: "
         << cylinder_filter_->id();
    fatal_error(mssg.str());
  }

  StaticVector3 cylinder_shape = cylinder_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), cylinder_shape.begin(),
                     cylinder_shape.end());
  axial_direction_ = cylinder_filter_->get_axial_direction();

  // last-dimension is sum of zernike and legendre orders to store the
  // coefficients of each polynomals. First will be zernike coefficients and
  // after will be the legendre coefficients.
  tally_shape.push_back((zr_order_ + 1) + (legen_order_ + 1));

  // reallocate and fill with zeros for the taly avg, gen-score and varaince
  tally_avg_.resize(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.resize(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.resize(tally_shape);
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

  if ((quantity_.type == Quantity::Type::Source ||
       quantity_.type == Quantity::Type::RealSource ||
       quantity_.type == Quantity::Type::ImagSource) &&
      estimator_ != Estimator::Source) {
    std::stringstream mssg;
    mssg << "Tally " << tally_name_
         << " has a source-like qantity but does not use a source estimator.";
    fatal_error(mssg.str());
  }

  StaticVector3 cylinder_shape = cylinder_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), cylinder_shape.begin(),
                     cylinder_shape.end());
  axial_direction_ = cylinder_filter_->get_axial_direction();

  // last-dimension is for the orders, since in this constructor, there is no
  // legendre evaluation. the length of the dimesion will be zr-order + 1
  tally_shape.push_back(zr_order_ + 1);

  // in-this particular constructor, there will be no legendre-order
  check_for_legendre = false;

  // reallocate and fill with zeros for the taly avg, gen-score and varaince
  tally_avg_.resize(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.resize(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.resize(tally_shape);
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
    // Not inside any cylinder bin. Don't score.
    return;
  }
  indices.insert(indices.end(), cylinder_index.begin(), cylinder_index.end());

  // add one dimension for the different orders of polynomials.
  // First loop over 0 to zernike-order for the zernike polynomial, then
  // loop over 0 + (zernike-order + 1) to legendre-order + (zernike-order + 1)
  const std::size_t FET_index = indices.size();
  indices.push_back(0);

  // Get the particle base score for the collision
  const double Et = mat.Et(p.E());
  const double collision_score =
      particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat) / Et;

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
    tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
  }

  // second- score for the legendre
  if (check_for_legendre == true) {
    const double zmin = cylinder_filter_->z_min(cylinder_index);
    const double inv_dz = cylinder_filter_->inv_dz();
    double z = tktr.r().z();
    if (axial_direction_ == CylinderFilter::Orientation::Y) {
      z = tktr.r().y();
    } else if (axial_direction_ == CylinderFilter::Orientation::X) {
      z = tktr.r().x();
    }
    const double scaled_z = 2. * (z - zmin) * inv_dz - 1.0;
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    for (std::size_t i = 0; i <= legen_order_; i++) {
      // score for i-th order's basis function
      if (i > 0) {
        // recursive relation to evaluate the legendre
        p2 = (scaled_z * static_cast<double>(2 * i - 1) * p1 -
              static_cast<double>(i - 1) * p0) /
             static_cast<double>(i);
        p0 = p1;
        p1 = p2;
      }
      beta_n = collision_score * p2;
      indices[FET_index] = i + zr_order_ + 1;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
    }
  }
}

void ZernikeFET::score_source(const BankedParticle& p) {
  const Position r = p.r;

  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_filter_) {
    std::optional<std::size_t> E_indx = energy_filter_->get_index(p.E);
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    indices.push_back(E_indx.value());
  }

  // get the positional indexes from cylinder_filter
  StaticVector3 cylinder_index = cylinder_filter_->get_position_index(r);
  if (cylinder_index.empty()) {
    // Not inside any cylinder bin. Don't score.
    return;
  }
  indices.insert(indices.end(), cylinder_index.begin(), cylinder_index.end());

  // Add one dimension for the different polynomial orders.
  // First loop over 0 to zernike-order for zernike polynomials.
  // After, loop over 0 + (zernike-order + 1) to
  // legendre-order + (zernike-order + 1) for the legendre polynomials.
  const std::size_t FET_index = indices.size();
  indices.push_back(0);

  // Get the particle base score for the source
  const double source_score = particle_base_score(p.E, p.wgt, p.wgt2, nullptr);

  // varaible for scoring
  double beta_n;

  // first- score for the zernike
  std::pair<double, double> scaled_r_and_theta =
      cylinder_filter_->get_scaled_radius_and_angle(cylinder_index, r);
  const double scaled_r = scaled_r_and_theta.first;
  const double theta = scaled_r_and_theta.second;
  const std::vector<double> zr_value =
      zr_polynomial_.evaluate_zernikes(scaled_r, theta);

  for (std::size_t i = 0; i <= zr_order_; i++) {
    // score for i-th order's basis function
    beta_n = source_score * zr_value[i];
    indices[FET_index] = i;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
  }

  // second- score for the legendre
  if (check_for_legendre == true) {
    const double zmin = cylinder_filter_->z_min(cylinder_index);
    const double inv_dz = cylinder_filter_->inv_dz();
    double z = r.z();
    if (axial_direction_ == CylinderFilter::Orientation::Y) {
      z = r.y();
    } else if (axial_direction_ == CylinderFilter::Orientation::X) {
      z = r.x();
    }
    const double scaled_z = 2. * (z - zmin) / inv_dz - 1.0;
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    for (std::size_t i = 0; i <= legen_order_; i++) {
      // score for i-th order's basis function
      if (i > 0) {
        // recursive relation to evaluate the legendre
        p2 = (scaled_z * static_cast<double>(2 * i - 1) * p1 -
              static_cast<double>(i - 1) * p0) /
             static_cast<double>(i);
        p0 = p1;
        p1 = p2;
      }
      beta_n = source_score * p2;
      indices[FET_index] = i + zr_order_ + 1;
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
    }
  }
}

double ZernikeFET::evaluate(const Position& r, const double& E) const {
  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_filter_) {
    std::optional<std::size_t> E_indx = energy_filter_->get_index(E);
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. return 0.;
      return 0.;
    }
    indices.push_back(E_indx.value());
  }

  // get the positional indexes from cylinder_filter
  StaticVector3 cylinder_index = cylinder_filter_->get_position_index(r);
  if (cylinder_index.empty()) {
    // Not inside any cylinder bin. return 0.
    return 0.;
  }
  indices.insert(indices.end(), cylinder_index.begin(), cylinder_index.end());

  // add one dimension for the different orders in polynomial
  // first loop over 0 to zernike-order for zernike
  // second loop over 0 + (zernike-order + 1) to legendre-order + (zernike-order
  // + 1)
  const std::size_t FET_index = indices.size();
  indices.push_back(0);

  double tally_value = 1.;
  // first evaluate the zernike
  // get the scaled-r and theta
  std::pair<double, double> scaled_r_and_theta =
      cylinder_filter_->get_scaled_radius_and_angle(cylinder_index, r);
  const double scaled_r = scaled_r_and_theta.first;
  const double theta = scaled_r_and_theta.second;
  // get the zernike-polynomial values upto given order
  const std::vector<double> zr_value =
      zr_polynomial_.evaluate_zernikes(scaled_r, theta);
  // loop over all the zernike-order
  double tally_zr = 0.;
  for (std::size_t order = 0; order <= zr_order_; order++) {
    indices[FET_index] = order;
    const double kn = zr_polynomial_.orthonormalization_constant(
        order);  // pi is not multiplied here
    tally_zr += kn * tally_avg_.element(indices.begin(), indices.end()) *
                zr_value[order];
  }
  tally_value *= tally_zr;

  // second- evaluate the legendre if exist
  if (check_for_legendre) {
    // get the inverse of cylinder length for scaled-z and volume calculation
    const double inv_dz_ = cylinder_filter_->inv_dz();
    // get the zmin and correct z for scaled_z
    const double zmin = cylinder_filter_->z_min(cylinder_index);
    double z = r.z();
    if (axial_direction_ == CylinderFilter::Orientation::Y) {
      z = r.y();
    } else if (axial_direction_ == CylinderFilter::Orientation::X) {
      z = r.x();
    }
    const double scaled_z = 2. * (z - zmin) * inv_dz_ - 1.0;

    // get the zero moment to normalise it
    indices[FET_index] = 0;
    tally_value /= tally_avg_.element(indices.begin(), indices.end());

    // loop over each legendre order for fet evaluation
    double tally_legndre = 0.;
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    for (std::size_t order = 0; order <= legen_order_; order++) {
      indices[FET_index] = order + zr_order_ + 1;
      // recursive relation to evaluate the legendre
      if (order > 0) {
        p2 = (scaled_z * static_cast<double>(2 * order - 1) * p1 -
              static_cast<double>(order - 1) * p0) /
             static_cast<double>(order);
        p0 = p1;
        p1 = p2;
      }
      tally_legndre += static_cast<double>(2 * order + 1) *
                       tally_avg_.element(indices.begin(), indices.end()) * p2;
    }
    tally_value *= tally_legndre;
  }
  return tally_value;
}

std::vector<double> ZernikeFET::evaluate(
    const std::vector<std::pair<Position, double>> r_E) const {
  // store the tally values
  std::vector<double> tallied_values;
  tallied_values.reserve(r_E.size());

  // loop over the positions to get first the energy and position index
  for (std::size_t i = 0; i < r_E.size(); i++) {
    StaticVector6 indices;
    const Position r = r_E[i].first;
    const double E = r_E[i].second;
    // get the energy-index, if energy-filter exists
    if (energy_filter_) {
      std::optional<std::size_t> E_indx = energy_filter_->get_index(E);
      if (E_indx.has_value() == false) {
        // Not inside any energy bin. push_back 0.;
        tallied_values.push_back(0.);
        continue;
      }
      indices.push_back(E_indx.value());
    }

    // get the positional indexes from cylinder_filter
    StaticVector3 cylinder_index = cylinder_filter_->get_position_index(r);
    if (cylinder_index.empty()) {
      // Not inside any cylinder bin. push_back 0.
      tallied_values.push_back(0.);
      continue;
    }
    indices.insert(indices.end(), cylinder_index.begin(), cylinder_index.end());

    // add one dimension for the different orders in polynomial
    // first loop over 0 to zernike-order for zernike
    // second loop over 0 + (zernike-order + 1) to legendre-order +
    // (zernike-order + 1)
    const std::size_t FET_index = indices.size();
    indices.push_back(0);

    double tally_value = 1.;
    // first evaluate the zernike
    // get the scaled-r and theta
    std::pair<double, double> scaled_r_and_theta =
        cylinder_filter_->get_scaled_radius_and_angle(cylinder_index, r);
    const double scaled_r = scaled_r_and_theta.first;
    const double theta = scaled_r_and_theta.second;
    // get the zernike-polynomial values upto given order
    const std::vector<double> zr_value =
        zr_polynomial_.evaluate_zernikes(scaled_r, theta);
    // loop over all the zernike-order
    double tally_zr = 0.;
    for (std::size_t order = 0; order <= zr_order_; order++) {
      indices[FET_index] = order;
      const double kn = zr_polynomial_.orthonormalization_constant(
          order);  // pi is not multiplied here
      tally_zr += kn * tally_avg_.element(indices.begin(), indices.end()) *
                  zr_value[order];
    }
    tally_value *= tally_zr;

    // get the inverse of cylinder length for scaled-z and volume calculation
    const double inv_dz_ = cylinder_filter_->inv_dz();
    // second- evaluate the legendre if exist
    if (check_for_legendre) {
      // get the zmin and correct z for scaled_z
      const double zmin = cylinder_filter_->z_min(cylinder_index);
      double z = r.z();
      if (axial_direction_ == CylinderFilter::Orientation::Y) {
        z = r.y();
      } else if (axial_direction_ == CylinderFilter::Orientation::X) {
        z = r.x();
      }
      const double scaled_z = 2. * (z - zmin) * inv_dz_ - 1.0;

      // get the zero moment to normalise it
      indices[FET_index] = 0;
      tally_value /= tally_avg_.element(indices.begin(), indices.end());

      // loop over each legendre order for fet evaluation
      double tally_legndre = 0.;
      double p0 = 1.;
      double p1 = 1.;
      double p2 = 1.;
      for (std::size_t order = 0; order <= legen_order_; order++) {
        indices[FET_index] = order + zr_order_ + 1;
        // recursive relation to evaluate the legendre
        if (order > 0) {
          p2 = (scaled_z * static_cast<double>(2 * order - 1) * p1 -
                static_cast<double>(order - 1) * p0) /
               static_cast<double>(order);
          p0 = p1;
          p1 = p2;
        }
        tally_legndre += static_cast<double>(2 * order + 1) *
                         tally_avg_.element(indices.begin(), indices.end()) *
                         p2;
      }
      tally_value *= tally_legndre;
    }
    tallied_values.push_back(tally_value);
  }
  return tallied_values;
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
  if (quantity_.type == Quantity::Type::MT) {
    tally_grp.createAttribute("mt", quantity_.mt);
  }

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
  this->var_to_std_on_mean();

  // Add data sets for the average and the standard deviation
  std::vector<std::size_t> shape(tally_avg_.shape().begin(),
                                 tally_avg_.shape().end());
  auto avg_dset = tally_grp.createDataSet<double>("avg", H5::DataSpace(shape));
  avg_dset.write_raw(tally_avg_.data());

  auto std_dset = tally_grp.createDataSet<double>("std", H5::DataSpace(shape));
  std_dset.write_raw(tally_var_.data());
}

// make the tally-zernike
std::shared_ptr<ZernikeFET> make_zernike_fet(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || node["name"].IsScalar() == false) {
    fatal_error("No valid name is provided on tally.");
  }
  std::string name = node["name"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    fatal_error("No quantity is given for tally " + name + ".");
  }
  given_quantity = node["quantity"].as<std::string>();
  Quantity quant = read_quantity(node, name);
  const bool source_like = (quant.type == Quantity::Type::Source ||
                            quant.type == Quantity::Type::RealSource ||
                            quant.type == Quantity::Type::ImagSource);

  std::string estimator_name = "collision";
  if (source_like) {
    estimator_name = "source";
  }

  // Check the estimator is given or not. We default to collision estimators.
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
  } else if (estimator_name == "source") {
    estimator = Estimator::Source;
  } else {
    std::stringstream mssg;
    mssg << "The tally " << name << " was provided with an unkown estimator \""
         << estimator_name << "\".";
    fatal_error(mssg.str());
  }

  if (source_like && estimator != Estimator::Source) {
    std::stringstream mssg;
    mssg << "Tally " << name
         << " has a source-like quantity but not a source estimator.";
    fatal_error(mssg.str());
  }

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
    fatal_error("A valid zernike-order was not given for " + name + " tally.");
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
    // check if the legendre exist, then make sure, it doesn't have a
    // infinte-cylinder
    if (cylinder_filter->is_infinite_cylinder() == true) {
      std::stringstream mssg;
      mssg << "Tally " << name
           << "cannot evaluate the legendre-fet over a infinte-cylinder with "
              "position-filter id: "
           << position_id;
      fatal_error(mssg.str());
    }
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