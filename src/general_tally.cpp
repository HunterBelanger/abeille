#include <tallies/general_tally.hpp>
#include <utils/error.hpp>

GeneralTally::GeneralTally(std::shared_ptr<PositionFilter> position_filter,
                           std::shared_ptr<EnergyFilter> energy_in,
                           Quantity quantity, Estimator estimator,
                           std::string name_)
    : ITally(quantity, estimator, name_),
      position_filter_(position_filter),
      energy_in_(energy_in) {
  // flag variable to predict weather position and energy filter exists or not
  bool check_ifany_filter = true;

  // get the shape or dimension of position filter
  StaticVector3 position_shape_;
  if (position_filter_) {
    position_shape_ = position_filter_->get_shape();
    check_ifany_filter = false;
  }

  StaticVector4 tally_dimensions_(position_shape_.begin(),
                                  position_shape_.end());

  // get the size of the enery_filter
  if (energy_in_) {
    size_t ne = energy_in_->size();
    tally_dimensions_.insert(tally_dimensions_.begin(), ne);

    check_ifany_filter = check_ifany_filter && false;
  }

  // if any of the filter is not present, then give fatal_error
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

  // get the indexes for the positions, if exist
  StaticVector3 position_indexes;
  if (position_filter_) {
    position_indexes = position_filter_->get_indices(
        tktr);  // it will provde the reduce dimensions
  }
  StaticVector4 indexes(position_indexes.begin(), position_indexes.end());

  // get the index for the energy if exist
  std::size_t index_E;
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value()) {
      index_E = E_indx.value();
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

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr,
                                double d_flight, MaterialHelper& mat) {
  if (!(estimator_ == GeneralTally::Estimator::TrackLength)) {
    return;
  }

  // get the energy index
  std::size_t index_E;
  if (energy_in_) {
    // get the energy_index if we in the bounds
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value()) {
      index_E = E_indx.value();
    } else {
      return;
    }
  }

  double flight_score = particle_base_score(p, mat);

  // get the all the bins' index along with the distance travelled by the
  // particle
  std::vector<TracklengthDistance> pos_indexes_ =
      position_filter_->get_indices_tracklength(trkr, d_flight);

  for (size_t iter = 0; iter < pos_indexes_.size(); iter++) {
    StaticVector4 all_indexes_(pos_indexes_[iter].indexes_.begin(),
                               pos_indexes_[iter].indexes_.end());
    const double d_ = pos_indexes_[iter].distance_in_bin;

    if (energy_in_) {
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
    // track-length can be used for flux
    if (!(given_quantity == "flux")) {
      fatal_error(
          "For tally \"" + general_tally_name + "\", the quantity \"" +
          given_quantity +
          "\" cannot be evaluated with the selected transport operator.");
    }

  itally_genral_tally = std::make_shared<GeneralTally>(
      position_filter_, energy_filter_, quant,
      GeneralTally::Estimator::TrackLength, general_tally_name);

  if (itally_genral_tally == nullptr)
    fatal_error("Incorrect \"estimator\" is given for " + general_tally_name +
                " tally.");

  return itally_genral_tally;
}
