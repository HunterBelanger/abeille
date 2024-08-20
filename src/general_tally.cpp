#include <tallies/general_tally.hpp>
#include <tallies/tallies.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector4 = boost::container::static_vector<size_t, 4>;

GeneralTally::GeneralTally(std::shared_ptr<PositionFilter> position_filter,
                           std::shared_ptr<EnergyFilter> energy_in,
                           Quantity quantity, Estimator estimator,
                           std::string name)
    : ITally(quantity, estimator, name),
      position_filter_(position_filter),
      energy_in_(energy_in) {
  StaticVector4 tally_shape;
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    // If no filters, make sure we have at least 1 element in the tally
    tally_shape.push_back(1);
  } else {
    // get the size of the enery_filter
    if (energy_in_) {
      size_t ne = energy_in_->size();
      tally_shape.push_back(ne);
    }

    // get the shape or dimension of position filter
    if (position_filter_) {
      StaticVector3 position_shape = position_filter_->get_shape();
      tally_shape.insert(tally_shape.end(), position_shape.begin(),
                         position_shape.end());
    }
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

  tally_avg_.resize(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.resize(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.resize(tally_shape);
  tally_var_.fill(0.0);
}

void GeneralTally::score_collision(const Particle& p, const Tracker& trkr,
                                   MaterialHelper& mat) {
  StaticVector4 indices;
  // if both energy and position filter don't exist, then index is 0
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    indices.push_back(0);
  } else {
    // get the index for the energy if exist
    if (energy_in_) {
      std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
      if (E_indx.has_value() == false) {
        // Not inside any energy bin. Don't score.
        return;
      }

      indices.push_back(E_indx.value());
    }

    // get the indices for the positions, if exist
    StaticVector3 position_indices;
    if (position_filter_) {
      // It will provde the reduce dimensions
      position_indices = position_filter_->get_indices(trkr);

      if (position_indices.empty()) {
        // Not inside any position bin. Don't score.
        return;
      }
    }
    indices.insert(indices.end(), position_indices.begin(),
                   position_indices.end());
  }

  const double Et = mat.Et(p.E());
  const double collision_score =
      particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat) / Et;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += collision_score;
}

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr,
                                double d_flight, MaterialHelper& mat) {
  // if position filter and energy-filter are not exit, then score and return
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    const double flight_score =
        particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat);
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_(0) += d_flight * flight_score;

    return;
  }

  // get the energy index
  std::size_t index_E;
  if (energy_in_) {
    // get the energy_index if we are in the bounds
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value() == false) {
      // Not inside any energy bin. Don't score.
      return;
    }

    index_E = E_indx.value();
  }

  const double flight_score =
      particle_base_score(p.E(), p.wgt(), p.wgt2(), &mat);

  if (position_filter_ == nullptr) {
    StaticVector4 all_indices;
    if (energy_in_) {
      all_indices.push_back(index_E);
    }

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_.element(all_indices.begin(), all_indices.end()) +=
        d_flight * flight_score;

    return;
  }

  // Get the all the bin indices along with the distance traveled by the
  // particle
  std::vector<TracklengthDistance> position_indices =
      position_filter_->get_indices_tracklength(trkr, d_flight);

  for (std::size_t iter = 0; iter < position_indices.size(); iter++) {
    StaticVector4 all_indices;
    if (energy_in_) {
      all_indices.push_back(index_E);
    }
    all_indices.insert(all_indices.end(), position_indices[iter].index.begin(),
                       position_indices[iter].index.end());

    const double d = position_indices[iter].distance;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_.element(all_indices.begin(), all_indices.end()) +=
        d * flight_score;
  }
}

void GeneralTally::score_source(const BankedParticle& p) {
  // Initialize a tracker
  Tracker trkr(p.r, p.u);

  StaticVector4 indices;
  // if both energy and position filter don't exist, then index is 0
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    indices.push_back(0);
  } else {
    // get the index for the energy if exist
    if (energy_in_) {
      std::optional<std::size_t> E_indx = energy_in_->get_index(p.E);
      if (E_indx.has_value() == false) {
        // Not inside any energy bin. Don't score.
        return;
      }

      indices.push_back(E_indx.value());
    }

    // get the indices for the positions, if exist
    StaticVector3 position_indices;
    if (position_filter_) {
      // It will provde the reduce dimensions
      position_indices = position_filter_->get_indices(trkr);

      if (position_indices.empty()) {
        // Not inside any position bin. Don't score.
        return;
      }
    }
    indices.insert(indices.end(), position_indices.begin(),
                   position_indices.end());
  }

  const double source_score = particle_base_score(p.E, p.wgt, p.wgt2, nullptr);

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_.element(indices.begin(), indices.end()) += source_score;
}

// method to evaluate the tally
double GeneralTally::evaluate(const Position& r, const double& E) const {
  StaticVector4 indices;
  // if both energy and position filter don't exist, then index is 0
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    indices.push_back(0);
  } else {
    // get the index for the energy if exist
    if (energy_in_) {
      std::optional<std::size_t> E_indx = energy_in_->get_index(E);
      if (E_indx.has_value() == false) {
        // Not inside any energy bin. Therefore, return 0.
        return 0.;
      }
      indices.push_back(E_indx.value());
    }
    // create a temporary tracker
    Direction u;
    Tracker trkr(r, u);
    // get the indices for the positions, if exist
    StaticVector3 position_indices;
    if (position_filter_) {
      // It will provde the reduce dimensions
      position_indices = position_filter_->get_indices(trkr);

      if (position_indices.empty()) {
        // Not inside any position bin. Therefore, return 0.
        return 0.;
      }
      indices.insert(indices.end(), position_indices.begin(),
                     position_indices.end());
    }
  }

  return tally_avg_.element(indices.begin(), indices.end());
}

std::vector<double> GeneralTally::evaluate(
    const std::vector<std::pair<Position, double>> r_E) const {
  std::vector<double> tallied_values;
  tallied_values.reserve(r_E.size());

  for (std::size_t i = 0; i <= r_E.size(); i++) {
    StaticVector4 indices;
    const Position r = r_E[i].first;
    const double E = r_E[i].second;
    // if both energy and position filter don't exist, then index is 0
    if (position_filter_ == nullptr && energy_in_ == nullptr) {
      indices.push_back(0);
    } else {
      // get the index for the energy if exist
      if (energy_in_) {
        std::optional<std::size_t> E_indx = energy_in_->get_index(E);
        if (E_indx.has_value() == false) {
          // Not inside any energy bin. push_back 0 and continue
          tallied_values.push_back(0.);
          continue;
        }
        indices.push_back(E_indx.value());
      }
      // create a temporary tracker
      Direction u;
      Tracker trkr(r, u);
      // get the indices for the positions, if exist
      StaticVector3 position_indices;
      if (position_filter_) {
        // It will provde the reduce dimensions
        position_indices = position_filter_->get_indices(trkr);

        if (position_indices.empty()) {
          // Not inside any position bin. Therefore, return 0.
          tallied_values.push_back(0.);
          continue;
        }
        indices.insert(indices.end(), position_indices.begin(),
                       position_indices.end());
      }
    }
    const double tally_value =
        tally_avg_.element(indices.begin(), indices.end());
    tallied_values.push_back(tally_value);
  }
  
  return tallied_values;
}

void GeneralTally::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + tally_name_);

  // Save the type
  tally_grp.createAttribute("type", "general-tally");

  // Save the quantity
  tally_grp.createAttribute("quantity", quantity_str());
  if (quantity_.type == Quantity::Type::MT) {
    tally_grp.createAttribute("mt", quantity_.mt);
  }

  // Save the estimator
  tally_grp.createAttribute("estimator", estimator_str());

  // Save energy-in filter id
  if (energy_in_) {
    tally_grp.createAttribute("energy-filter", energy_in_->id());
  }

  // Save position filter id
  if (position_filter_) {
    tally_grp.createAttribute("position-filter", position_filter_->id());
  }

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

std::shared_ptr<GeneralTally> make_general_tally(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || node["name"].IsScalar() == false) {
    fatal_error("No valid name is provided on tally.");
  }
  std::string name = node["name"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    fatal_error("No valid quantity entry is given for tally " + name + ".");
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
  } else if (estimator_name == "source") {
    estimator = Estimator::Source;
  } else {
    fatal_error("Invalid estimator is given for tally " + name + ".");
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

  // Get the position filter
  std::shared_ptr<PositionFilter> position_filter = nullptr;
  if (node["position-filter"] && node["position-filter"].IsScalar()) {
    std::size_t position_id = node["position-filter"].as<std::size_t>();
    position_filter = tallies.get_position_filter(position_id);
    if (position_filter == nullptr) {
      std::stringstream mssg;
      mssg << "For tally " << name << ", cannot find position filter with id "
           << position_id << ".";
      fatal_error(mssg.str());
    }
  }

  // For the general tally
  std::shared_ptr<GeneralTally> tally = std::make_shared<GeneralTally>(
      position_filter, energy_filter, quant, estimator, name);

  return tally;
}
