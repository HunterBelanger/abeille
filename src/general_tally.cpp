#include <tallies/general_tally.hpp>
#include <tallies/itally.hpp>
#include <utils/error.hpp>

#include <boost/container/static_vector.hpp>

using StaticVector4 = boost::container::static_vector<size_t, 4>;

GeneralTally::GeneralTally(std::shared_ptr<PositionFilter> position_filter,
                           std::shared_ptr<EnergyFilter> energy_in,
                           Quantity quantity, Estimator estimator,
                           std::string name_)
    : ITally(quantity, estimator, name_),
      position_filter_(position_filter),
      energy_in_(energy_in) {
        
  bool check_ifany_filter = true;

  StaticVector3 position_shape_;
  if (position_filter_) {
    position_shape_ = position_filter_->get_dimension();
    check_ifany_filter = false;
  }

  StaticVector4 tally_dimensions_(position_shape_.begin(), position_shape_.end());

  if (energy_in_) {
    size_t ne = energy_in_->size();
    tally_dimensions_.insert(tally_dimensions_.begin(), ne);

    check_ifany_filter = check_ifany_filter && false;
  }

  if (check_ifany_filter == true) {
    fatal_error("for the " + name_ + " tally, no filter is provided.");
  }

  tally_avg.reallocate(tally_dimensions_);
  tally_avg.fill(0.0);

  tally_gen_score.reallocate(tally_dimensions_);
  tally_gen_score.fill(0.0);

  tally_var.reallocate(tally_dimensions_);
  tally_var.fill(0.0);
}

void GeneralTally::score_collision(const Particle& p, const Tracker& tktr,
                                   MaterialHelper& mat) {
  if (!(estimator_ == GeneralTally::Estimator::Collision)) {
    return;
  }

  std::size_t index_E;
  StaticVector3 position_indexes;

  if (position_filter_) {
    position_indexes = position_filter_->get_indices(
        tktr);  // it will provde the reduce dimensions
  }
  
  StaticVector4 indexes(position_indexes.begin(), position_indexes.end());

  if (energy_in_) {
    if (energy_in_->get_index(p.E(), index_E)) {
      indexes.insert(indexes.begin(), index_E);

    } else {
      return;
    }
  }

  if (indexes.empty()) {
    return;
  }

  const double Et = mat.Et(p.E());
  const double collision_score = particle_base_score(p, mat) / Et;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score(indexes) += collision_score;
}

// void GeneralTally::score_flight(const Particle& p,
//                                 double d_flight ,MaterialHelper& mat ){

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr,
                                double d_flight, MaterialHelper& mat) {
  if (!(estimator_ == GeneralTally::Estimator::TrackLength)) {
    return;
  }

  std::size_t index_E;
  bool any_energy_in_ = false;  // Should be true when the energy_in_ is there

  if (energy_in_) {
    // get the energy_index if we in the bounds
    if (energy_in_->get_index(p.E(), index_E)) {
      // for track-length, it will done with help of bool
      any_energy_in_ = true;
    } else {
      return;
    }
  }

  double flight_score = particle_base_score(p, mat);

  std::vector<TracklengthDistance> pos_indexes_ =
      position_filter_->get_indices_tracklength(trkr, d_flight);

  for (size_t iter = 0; iter < pos_indexes_.size(); iter++) {
    std::vector<size_t> all_indexes_ = pos_indexes_[iter].indexes_;
    const double d_ = pos_indexes_[iter].distance_in_bin;

    if (any_energy_in_ == false) {
      all_indexes_.insert(all_indexes_.begin(), index_E);
    }

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score(all_indexes_) += d_ * flight_score;
  }
}

std::shared_ptr<GeneralTally> make_general_tally(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"]) {
    fatal_error("No valid name is provided.");
  }
  std::string general_tally_name = node["name"].as<std::string>();

  // Check the estimator is given or not.
  if (!node["estimator"]) {
    fatal_error(
        "No, estimator is given for \"" + general_tally_name +
        "\", therefore, collision estimator will be taken, by default.");
  }
  std::string estimator_name = node["estimator"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"]) {
    fatal_error("No quantity is given for " + general_tally_name + " tally.");
  }
  given_quantity = node["quantity"].as<std::string>();
  GeneralTally::Quantity quant;

  bool found_quantity = false;
  if (given_quantity == "flux") {
    quant = GeneralTally::Quantity::Flux;
    found_quantity = true;
  }

  if (given_quantity == "fission") {
    quant = GeneralTally::Quantity::Fission;
    found_quantity = true;
  }
  if (given_quantity == "absorption") {
    quant = GeneralTally::Quantity::Absorption;
    found_quantity = true;
  }
  if (given_quantity == "elastic") {
    quant = GeneralTally::Quantity::Elastic;
    found_quantity = true;
  }

  if (found_quantity == false) {
    fatal_error("For " + general_tally_name +
                " tally, a unknown quantity is given.");
  }

  // Get the enrgy bounds, if any is given
  std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;
  if (node["energy-bounds"]) {
    std::vector<double> energy_bounds =
        node["energy-bounds"].as<std::vector<double>>();
    energy_filter_ = std::make_shared<EnergyFilter>(energy_bounds);
  }

  // Get the position filter
  std::shared_ptr<PositionFilter> position_filter_ = make_position_filter(node);

  // For the general tally
  std::shared_ptr<GeneralTally> itally_genral_tally = nullptr;
  if (estimator_name == "collision") {
    itally_genral_tally = std::make_shared<GeneralTally>(
        position_filter_, energy_filter_, quant,
        GeneralTally::Estimator::Collision, general_tally_name);
  }

  if (estimator_name == "track-length")
    itally_genral_tally = std::make_shared<GeneralTally>(
        position_filter_, energy_filter_, quant,
        GeneralTally::Estimator::TrackLength, general_tally_name);

  if (itally_genral_tally == nullptr)
    fatal_error("Incorrect \"estimator\" is given for " + general_tally_name +
                " tally.");

  return itally_genral_tally;
}
