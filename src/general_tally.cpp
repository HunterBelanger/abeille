#include <tallies/general_tally.hpp>
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
  // at least one filter must exists.
  if ( position_filter_ == nullptr && energy_in_ == nullptr){
      fatal_error("for the " + name_ + " tally, no filter is provided.");
  }

  // get the shape or dimension of position filter
  StaticVector3 position_shape_;
  if (position_filter_) {
    position_shape_ = position_filter_->get_shape();
  }

  StaticVector4 tally_dimensions_(position_shape_.begin(),
                                  position_shape_.end());

  // get the size of the enery_filter
  if (energy_in_) {
    size_t ne = energy_in_->size();
    tally_dimensions_.insert(tally_dimensions_.begin(), ne);
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
  if (estimator_ != GeneralTally::Estimator::Collision) {
    return;
  }

  // get the indices for the positions, if exist
  StaticVector3 position_indices;
  if (position_filter_) {
    position_indices = position_filter_->get_indices(
        tktr);  // it will provde the reduce dimensions
  }
  StaticVector4 indices(position_indices.begin(), position_indices.end());

  // get the index for the energy if exist
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    if (E_indx.has_value()) {
      std::size_t index_E = E_indx.value();
      indices.insert( indices.begin(), index_E );
    } else {
      // Not inside any energy bin. Don't score.
      return;
    }
  }

  if (indices.empty()) {
    return;
  }

  const double Et = mat.Et(p.E());
  const double collision_score = particle_base_score(p, mat) / Et;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  tally_gen_score(indices) += collision_score;
}

void GeneralTally::score_flight(const Particle& p, const Tracker& trkr,
                                double d_flight, MaterialHelper& mat) {
  if (estimator_ != GeneralTally::Estimator::TrackLength) {
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
      // Not inside any energy bin. Don't score.
      return;
    }
  }

  double flight_score = particle_base_score(p, mat);

  // get the all the bins' index along with the distance travelled by the
  // particle
  std::vector<TracklengthDistance> pos_indices_ =
      position_filter_->get_indices_tracklength(trkr, d_flight);

  for (size_t iter = 0; iter < pos_indices_.size(); iter++) {
    StaticVector4 all_indices_(pos_indices_[iter].index.begin(),
                               pos_indices_[iter].index.end());
    
    if (energy_in_) {
      all_indices_.insert(all_indices_.begin(), index_E);
    }

    const double d_ = pos_indices_[iter].distance;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen_score(all_indices_) += d_ * flight_score;
  }
}

std::shared_ptr<GeneralTally> make_general_tally(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] || !node["name"].IsScalar()) {
    fatal_error("No valid name is provided.");
  }
  std::string general_tally_name = node["name"].as<std::string>();

  // Check the estimator is given or not.
  if (!node["estimator"] || !node["estimator"].IsScalar()) {
    fatal_error("Estimator is not given for \"" + general_tally_name + ".");
  }
  std::string estimator_name = node["estimator"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    fatal_error("No quantity is given for " + general_tally_name + " tally.");
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
    fatal_error("For " + general_tally_name +
                " tally, a unknown quantity is given.");
  }

  // Get the enrgy bounds, if any is given
  std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;
  if (node["energy-bounds"]) {
    if ( !node["energy-bounds"].IsScalar() ){
      fatal_error("Invalid energy-filter is given.");
    }
    energy_filter_ = make_energy_filter(node);
  }

  // Get the position filter
  std::shared_ptr<PositionFilter> position_filter_ = nullptr;
  if (node["position-filter"]) {
    if ( !node["position-filter"].IsScalar() ){
      fatal_error("Invalid position-filter is given.");
    }
    position_filter_ = make_position_filter(node);
  }

  // For the general tally
  std::shared_ptr<GeneralTally> itally_general_tally = nullptr;
  if (estimator_name == "collision") {
    itally_general_tally = std::make_shared<GeneralTally>(
        position_filter_, energy_filter_, quant,
        GeneralTally::Estimator::Collision, general_tally_name);
  } else if (estimator_name == "track-length") {
    itally_general_tally = std::make_shared<GeneralTally>(
        position_filter_, energy_filter_, quant,
        GeneralTally::Estimator::TrackLength, general_tally_name);

  } else {
    fatal_error("Invalid estimator is given for \"" + general_tally_name + ".");
  }

  return itally_general_tally;
}
