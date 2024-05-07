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

  tally_avg_.reallocate(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.reallocate(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.reallocate(tally_shape);
  tally_var_.fill(0.0);
}

void GeneralTally::score_collision(const Particle& p, const Tracker& tktr,
                                   MaterialHelper& mat) {
  StaticVector4 indices;
  // if both energy and position filter don't exist, then index is 0
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    indices.push_back(0);
  } else {
    // get the index for the energy if exist
    if (energy_in_) {
      std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
      if (E_indx.has_value()) {
        std::size_t index_E = E_indx.value();
        indices.push_back(index_E);
      } else {
        // Not inside any energy bin. Don't score.
        return;
      }
    }
    // get the indices for the positions, if exist
    StaticVector3 position_indices;
    if (position_filter_) {
      position_indices = position_filter_->get_indices(
          tktr);  // it will provde the reduce dimensions

      if (position_indices.empty()) {
        // Not inside any energy bin. Don't score.
        return;
      }
    }
    indices.insert(indices.end(), position_indices.begin(),
                   position_indices.end());
  }

  const double Et = mat.Et(p.E());
  const double collision_score = particle_base_score(p, mat) / Et;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score_(indices) += collision_score;
}

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr,
                                double d_flight, MaterialHelper& mat) {
  // if position filter and energy-filter are not exit, then score and return
  if (position_filter_ == nullptr && energy_in_ == nullptr) {
    const double flight_score = particle_base_score(p, mat);
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score_({0}) += d_flight * flight_score;

    return;
  }

  // get the energy index
  std::size_t index_E;
  if (energy_in_) {
    // get the energy_index if we are in the bounds
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value()) {
      index_E = E_indx.value();
    } else {
      // Not inside any energy bin. Don't score.
      return;
    }
  }

  const double flight_score = particle_base_score(p, mat);

  // get the all the bins' index along with the distance traveled by the
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
    tally_gen_score_(all_indices) += d * flight_score;
  }
}

void GeneralTally::write_tally() {
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

std::shared_ptr<GeneralTally> make_general_tally(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || !node["name"].IsScalar()) {
    fatal_error("No valid name is provided.");
  }
  std::string name = node["name"].as<std::string>();

  // Check the estimator is given or not. We default to collision estimators.
  std::string estimator_name = "collision";
  if (node["estimator"] && node["estimator"].IsScalar()) {
    estimator_name = node["estimator"].as<std::string>();
  } else if (node["estimator"]) {
    fatal_error("Invalid estimator entry is given for \"" + name + "\".");
  }

  GeneralTally::Estimator estimator;
  if (estimator_name == "collision") {
    estimator = GeneralTally::Estimator::Collision;
  } else if (estimator_name == "track-length") {
    estimator = GeneralTally::Estimator::TrackLength;
  } else {
    fatal_error("Invalid estimator is given for \"" + name + ".");
  }

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    fatal_error("No quantity is given for " + name + " tally.");
  }
  given_quantity = node["quantity"].as<std::string>();
  GeneralTally::Quantity quant;

  if (given_quantity == "flux") {
    quant = GeneralTally::Quantity::Flux;
  } else if (given_quantity == "fission") {
    quant = GeneralTally::Quantity::Fission;
  } else if (given_quantity == "absorption") {
    quant = GeneralTally::Quantity::Absorption;
  } else if (given_quantity == "elastic") {
    quant = GeneralTally::Quantity::Elastic;
  } else {
    fatal_error("For tally " + name + " unknown quantity \"" + given_quantity +
                "\" is given.");
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
