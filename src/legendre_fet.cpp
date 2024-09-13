#include <tallies/itally.hpp>
#include <tallies/legendre_fet.hpp>
#include <tallies/tallies.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<size_t, 6>;

LegendreFET::LegendreFET(std::shared_ptr<CartesianFilter> position_filter,
                         std::shared_ptr<EnergyFilter> energy_in,
                         std::vector<LegendreFET::Axis> axes,
                         std::vector<std::size_t> fet_order, Quantity quantity,
                         Estimator estimator, std::string name)
    : ITally(quantity, estimator, name),
      cartesian_filter_(position_filter),
      energy_in_(energy_in),
      axes_(axes.begin(), axes.end()),
      fet_order_(fet_order) {
  StaticVector6 tally_shape;
  // add the dimension for energy_in_ only if exist
  if (energy_in_) {
    std::size_t ne = energy_in_->size();
    tally_shape.push_back(ne);
  }

  // get the shape or dimensions for the cartesian_filter_
  // the legendre-fet requires the co-ordinate information

  if (cartesian_filter_ == nullptr) {
    fatal_error("LegendreFET has nullptr cartesian-filter.");
  }
  StaticVector3 position_shape = cartesian_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), position_shape.begin(),
                     position_shape.end());

  // for legendre, check the size of axis vector should be between than 3.
  if (axes_.size() <= 1 && axes_.size() >= 3)
    fatal_error("LegendreFET length of axes must be in [1,3].");

  if ((quantity_.type == Quantity::Type::Source ||
       quantity_.type == Quantity::Type::RealSource ||
       quantity_.type == Quantity::Type::ImagSource) &&
      estimator_ != Estimator::Source) {
    std::stringstream mssg;
    mssg << "Tally " << tally_name_
         << " has a source-like qantity but does not use a source estimator.";
    fatal_error(mssg.str());
  }

  // Add the dimension for the order of FET only
  // the dimension size will be sum of all order+1 for any given axis
  std::size_t fet_order_dimension = 0;
  for (std::size_t it_axis = 0; it_axis < fet_order_.size(); it_axis++) {
    fet_order_dimension += fet_order_[it_axis] + 1;
  }
  tally_shape.push_back(fet_order_dimension);

  // reallocate and fill with zeros for the tally avg, gen-score and variance
  tally_avg_.resize(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.resize(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.resize(tally_shape);
  tally_var_.fill(0.0);
}

void LegendreFET::score_collision(const Particle& p, const Tracker& tktr,
                                  MaterialHelper& mat) {
  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    indices.push_back(E_indx.value());
  }

  // get the cartisian_filter indices
  StaticVector3 position_index = cartesian_filter_->get_indices(tktr);
  if (position_index.empty()) {
    // No bin is found, don't score.
    return;
  }

  indices.insert(indices.end(), position_index.begin(), position_index.end());

  const double Et = mat.Et(p.E());

  const double collision_score =
      particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat) / Et;

  // add one dimesnion for FET-index
  // this dimension will loop over all the order for differnt axes
  const size_t FET_index = indices.size();
  indices.push_back(0);
  // index to iterate over the differnt order in axes
  std::size_t it_coeff = 0;

  // Variables for scoring
  double beta_n, scaled_loc;
  // Loop over the different axis and indexing is done
  for (std::size_t it_axis = 0; it_axis < axes_.size(); it_axis++) {
    // get the scaled x, y, or z for legendre polynomial
    switch (axes_[it_axis]) {
      case LegendreFET::Axis::X: {
        const double xmin_ = cartesian_filter_->x_min(position_index);
        const double inv_dx_ = cartesian_filter_->inv_dx(position_index);
        scaled_loc = 2. * (tktr.r().x() - xmin_) * inv_dx_ - 1.;
      } break;

      case LegendreFET::Axis::Y: {
        const double ymin_ = cartesian_filter_->y_min(position_index);
        const double inv_dy_ = cartesian_filter_->inv_dy(position_index);
        scaled_loc = 2. * (tktr.r().y() - ymin_) * inv_dy_ - 1.;
      } break;

      case LegendreFET::Axis::Z: {
        const double zmin_ = cartesian_filter_->z_min(position_index);
        const double inv_dz_ = cartesian_filter_->inv_dz(position_index);
        scaled_loc = 2. * (tktr.r().z() - zmin_) * inv_dz_ - 1.;
      }
    }

    // loop over differnt FET order
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    for (std::size_t i = 0; i <= fet_order_[it_axis]; i++) {
      if (i > 0) {
        // recursive relation to evaluate the legendre
        p2 = (scaled_loc * static_cast<double>(2 * i - 1) * p1 -
              static_cast<double>(i - 1) * p0) /
             static_cast<double>(i);
        p0 = p1;
        p1 = p2;
      }
      beta_n = collision_score * p2;
      indices[FET_index] = it_coeff;
      it_coeff++;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
    }
  }
}

void LegendreFET::score_source(const BankedParticle& p) {
  const Position r = p.r;

  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E);
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    indices.push_back(E_indx.value());
  }

  // get the cartisian_filter indices
  StaticVector3 position_index = cartesian_filter_->get_position_index(r);
  if (position_index.empty()) {
    // No bin is found, don't score.
    return;
  }

  indices.insert(indices.end(), position_index.begin(), position_index.end());

  const double source_score = particle_base_score(p.E, p.wgt, p.wgt2, nullptr);

  // add one dimesnion for FET-index
  // this dimension will loop over all the order for differnt axes
  const size_t FET_index = indices.size();
  indices.push_back(0);
  // index to iterate over the differnt order in axes
  std::size_t it_coeff = 0;

  // Variables for scoring
  double beta_n, scaled_loc;
  // Loop over the different axis and indexing is done
  for (std::size_t it_axis = 0; it_axis < axes_.size(); it_axis++) {
    // get the scaled x, y, or z for legendre polynomial
    switch (axes_[it_axis]) {
      case LegendreFET::Axis::X: {
        const double xmin_ = cartesian_filter_->x_min(position_index);
        const double inv_dx_ = cartesian_filter_->inv_dx(position_index);
        scaled_loc = 2. * (r.x() - xmin_) * inv_dx_ - 1.;
      } break;

      case LegendreFET::Axis::Y: {
        const double ymin_ = cartesian_filter_->y_min(position_index);
        const double inv_dy_ = cartesian_filter_->inv_dy(position_index);
        scaled_loc = 2. * (r.y() - ymin_) * inv_dy_ - 1.;
      } break;

      case LegendreFET::Axis::Z: {
        const double zmin_ = cartesian_filter_->z_min(position_index);
        const double inv_dz_ = cartesian_filter_->inv_dz(position_index);
        scaled_loc = 2. * (r.z() - zmin_) * inv_dz_ - 1.;
      }
    }

    // loop over differnt FET order
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    for (std::size_t i = 0; i <= fet_order_[it_axis]; i++) {
      // recursive relation to evaluate the legendre
      if (i > 0) {
        p2 = (scaled_loc * static_cast<double>(2 * i - 1) * p1 -
              static_cast<double>(i - 1) * p0) /
             static_cast<double>(i);
        p0 = p1;
        p1 = p2;
      }
      beta_n = source_score * p2;
      indices[FET_index] = it_coeff;
      it_coeff++;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_.element(indices.begin(), indices.end()) += beta_n;
    }
  }
}

double LegendreFET::evaluate(const Position& r, const double& E) const {
  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(E);
    if (E_indx.has_value() == false) {
      // Not inside any energy bin.return 0.
      return 0.;
    }
    indices.push_back(E_indx.value());
  }
  // get the cartisian_filter indices
  StaticVector3 position_index = cartesian_filter_->get_position_index(r);
  if (position_index.empty()) {
    // No bin is found, return 0.
    return 0.;
  }
  indices.insert(indices.end(), position_index.begin(), position_index.end());

  // add one dimesnion for FET-index
  // this dimension will loop over all the order for differnt axes
  const size_t FET_index = indices.size();
  indices.push_back(0);
  // index to iterate over the differnt order in axes
  std::size_t it_coeff = 0;

  double tally_value = 1.;  // for evaluation of tally at a given Position
  double scaled_loc;
  for (std::size_t it_axis = 0; it_axis < axes_.size(); it_axis++) {
    // check every existing axis and calculate for each and every order
    switch (axes_[it_axis]) {
      case LegendreFET::Axis::X: {
        // get the x-min of the position bin and calculate the scaled-x
        const double xmin_ = cartesian_filter_->x_min(position_index);
        // get the inverse of dx
        const double inv_dx_ = cartesian_filter_->inv_dx(position_index);
        scaled_loc = 2. * (r.x() - xmin_) * inv_dx_ - 1.;
      } break;

      case LegendreFET::Axis::Y: {
        // get the y-min of the position bin and calculate the scaled-y
        const double ymin_ = cartesian_filter_->y_min(position_index);
        // get the inverse of dy
        const double inv_dy_ = cartesian_filter_->inv_dy(position_index);
        scaled_loc = 2. * (r.y() - ymin_) * inv_dy_ - 1.;
      } break;

      case LegendreFET::Axis::Z: {
        // get the z-min of the position bin and calculate the scaled-z
        const double zmin_ = cartesian_filter_->z_min(position_index);
        // get the inverse of dz
        const double inv_dz_ = cartesian_filter_->inv_dz(position_index);
        scaled_loc = 2. * (r.z() - zmin_) * inv_dz_ - 1.;
      }
    }

    double fet_value = 0.;
    double p0 = 1.;
    double p1 = 1.;
    double p2 = 1.;
    // loop over each order and calculate the fet for respective axis
    for (std::size_t order = 0; order <= fet_order_[it_axis]; order++) {
      indices[FET_index] = it_coeff;
      it_coeff++;

      if (order > 0) {
        // recursive relation to evaluate the legendre
        p2 = (scaled_loc * static_cast<double>(2 * order - 1) * p1 -
              static_cast<double>(order - 1) * p0) /
             static_cast<double>(order);
        p0 = p1;
        p1 = p2;
      }
      fet_value += (2. * static_cast<double>(order) + 1.) *
                   tally_avg_.element(indices.begin(), indices.end()) * p2;
    }
    tally_value *= fet_value;

    // The normalisation is required and done by dividing zero-moment
    if (it_axis > 0) {
      indices[FET_index] = 0;
      tally_value /= tally_avg_.element(indices.begin(), indices.end());
    }
  }

  return tally_value;
}

std::vector<double> LegendreFET::evaluate(
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
    if (energy_in_) {
      std::optional<std::size_t> E_indx = energy_in_->get_index(E);
      if (E_indx.has_value() == false) {
        // Not inside any energy bin. push_back 0. in talied_values and
        // continue.
        tallied_values.push_back(0.);
        continue;
      }
      indices.push_back(E_indx.value());
    }
    // get the cartisian_filter indices
    StaticVector3 position_index = cartesian_filter_->get_position_index(r);
    if (position_index.empty()) {
      // No bin is found, push_back 0.
      tallied_values.push_back(0.);
      continue;
    }
    indices.insert(indices.end(), position_index.begin(), position_index.end());

    // add one dimesnion for FET-index
    // this dimension will loop over all the order for differnt axes
    const size_t FET_index = indices.size();
    indices.push_back(0);
    // index to iterate over the differnt order in axes
    std::size_t it_coeff = 0;

    double tally_value = 1.;  // for evaluation of tally at a given Position
    double scaled_loc;
    for (std::size_t it_axis = 0; it_axis < axes_.size(); it_axis++) {
      // check every existing axis and calculate for each and every order
      switch (axes_[it_axis]) {
        case LegendreFET::Axis::X: {
          // get the x-min of the position bin and calculate the scaled-x
          const double xmin_ = cartesian_filter_->x_min(position_index);
          // get the inverse of dx
          const double inv_dx_ = cartesian_filter_->inv_dx(position_index);
          scaled_loc = 2. * (r.x() - xmin_) * inv_dx_ - 1.;
        } break;

        case LegendreFET::Axis::Y: {
          // get the y-min of the position bin and calculate the scaled-y
          const double ymin_ = cartesian_filter_->y_min(position_index);
          // get the inverse of dy
          const double inv_dy_ = cartesian_filter_->inv_dy(position_index);
          scaled_loc = 2. * (r.y() - ymin_) * inv_dy_ - 1.;
        } break;

        case LegendreFET::Axis::Z: {
          // get the z-min of the position bin and calculate the scaled-z
          const double zmin_ = cartesian_filter_->z_min(position_index);
          // get the inverse of dz
          const double inv_dz_ = cartesian_filter_->inv_dz(position_index);
          scaled_loc = 2. * (r.z() - zmin_) * inv_dz_ - 1.;
        }
      }

      double fet_value = 0.;
      double p0 = 1.;
      double p1 = 1.;
      double p2 = 1.;
      // loop over each order and calculate the fet for respective axis
      for (std::size_t order = 0; order <= fet_order_[it_axis]; order++) {
        indices[FET_index] = it_coeff;
        it_coeff++;

        if (order > 0) {
          // recursive relation to evaluate the legendre
          p2 = (scaled_loc * static_cast<double>(2 * order - 1) * p1 -
                static_cast<double>(order - 1) * p0) /
               static_cast<double>(order);
          p0 = p1;
          p1 = p2;
        }
        fet_value += (2. * static_cast<double>(order) + 1.) *
                     tally_avg_.element(indices.begin(), indices.end()) * p2;
      }
      tally_value *= fet_value;

      // The normalisation is required and done by dividing zero-moment
      if (it_axis > 0) {
        indices[FET_index] = 0;
        tally_value /= tally_avg_.element(indices.begin(), indices.end());
      }
    }
    tallied_values.push_back(tally_value);
  }

  return tallied_values;
}

void LegendreFET::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + tally_name_);

  // Save the type
  tally_grp.createAttribute("type", "legendre-fet");

  // Save the energy-in filter
  if (energy_in_) {
    tally_grp.createAttribute("energy-filter", energy_in_->id());
  }

  // Save position filter
  tally_grp.createAttribute("position-filter", cartesian_filter_->id());

  // Save FET order
  tally_grp.createAttribute("order", fet_order_);

  // Save the axes
  std::vector<std::string> axes;
  for (const auto& ax : axes_) {
    switch (ax) {
      case Axis::X:
        axes.push_back("x");
        break;

      case Axis::Y:
        axes.push_back("y");
        break;

      case Axis::Z:
        axes.push_back("z");
        break;
    }
  }
  tally_grp.createAttribute("axes", axes);

  // Save the quantity
  tally_grp.createAttribute("quantity", quantity_str());
  if (quantity_.type == Quantity::Type::MT) {
    tally_grp.createAttribute("mt", quantity_.mt);
  }

  // Save the estimator
  tally_grp.createAttribute("estimator", estimator_str());

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

std::shared_ptr<LegendreFET> make_legendre_fet(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] && node["name"].IsScalar() == false) {
    fatal_error("No valid name is provided on tally.");
  }
  std::string name = node["name"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    fatal_error("The tally " + name +
                " has an invalid/nonexistent quantity entry.");
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
    fatal_error("Invalid estimator entry is given on tally " + name + ".");
  }
  Estimator estimator;
  if (estimator_name == "collision") {
    estimator = Estimator::Collision;
  } else if (estimator_name == "track-length") {
    estimator = Estimator::TrackLength;
    fatal_error(
        "On tally " + name +
        ", track-length estimator is not yet supported on legendre-fet.");
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

  // get the tallies instance
  auto& tallies = Tallies::instance();

  // Get the energy bounds, if any is given
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

  // Get the cartesian type position filter
  std::shared_ptr<CartesianFilter> cartesian_filter = nullptr;
  if (!node["position-filter"] || node["position-filter"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Legendre-FET tally " << name
         << " has invalid/nonexistent position-filter entry.";
    fatal_error(mssg.str());
  }
  std::size_t position_id = node["position-filter"].as<std::size_t>();
  cartesian_filter = tallies.get_cartesian_filter(position_id);
  if (cartesian_filter == nullptr) {
    std::stringstream mssg;
    mssg << "Legenre-FET tally " << name << " was provided position-filter id "
         << position_id << ". No cartesian filter with this id was found.";
    fatal_error(mssg.str());
  }

  // Get the legendre-fet order
  if (!node["order"] || !node["order"].IsSequence()) {
    std::stringstream mssg;
    mssg << "Legendre-FET tally " << name
         << " was not provided a valid order entry.";
    fatal_error(mssg.str());
  }
  std::vector<std::size_t> order = node["order"].as<std::vector<std::size_t>>();

  // Get the legendre-fet axes
  if (!node["axes"] || node["axes"].IsSequence() == false) {
    std::stringstream mssg;
    mssg << "Tally " << name << " does not have a valid axes list entry.";
    fatal_error(mssg.str());
  }
  std::vector<std::string> axis_strs =
      node["axes"].as<std::vector<std::string>>();

  std::vector<LegendreFET::Axis> axes;
  axes.reserve(axis_strs.size());
  for (const auto& c : axis_strs) {
    if (c == "X" || c == "x")
      axes.push_back(LegendreFET::Axis::X);
    else if (c == "Y" || c == "y")
      axes.push_back(LegendreFET::Axis::Y);
    else if (c == "Z" || c == "z")
      axes.push_back(LegendreFET::Axis::Z);
    else {
      std::stringstream mssg;
      mssg << "Tally " << name << " was provided invalid axis \"" << c << "\".";
      fatal_error(mssg.str());
    }
  }
  if (axes.empty() || axes.size() > 3) {
    std::stringstream mssg;
    mssg << "Tally " << name << " must have 1 to 3 axis entries.";
    fatal_error(mssg.str());
  }

  // Make the Legendre FET tally class
  std::shared_ptr<LegendreFET> tally = std::make_shared<LegendreFET>(
      cartesian_filter, energy_filter, axes, order, quant, estimator, name);

  return tally;
}