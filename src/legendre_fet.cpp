#include <tallies/legendre_fet.hpp>
#include <utils/legendre.hpp>
#include <tallies/itally.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<size_t, 6>;

LegendreFET::LegendreFET(std::shared_ptr<CartesianFilter> position_filter_,
                         std::shared_ptr<EnergyFilter> energy_in,
                         std::vector<LegendreFET::Axis> axes_,
                         size_t FET_order_, LegendreFET::Quantity quantity_,
                         LegendreFET::Estimator estimator_, std::string name_)
    : ITally(quantity_, estimator_, name_),
      cartesian_filter_(position_filter_),
      energy_in_(energy_in),
      axes(),
      FET_order(FET_order_) {
  // for legendre, check the size of axis vector should be between than 3.
  for ( auto& c : axes_){
    axes.push_back(c);
  }
  if (axes.size() <= 1 && axes.size() >= 3)
    fatal_error(
        "The no. of given axis for legendre-FET is not between 1 to 3 for "
        "tally " +
        name_ + ".");



  bool check_filter = true;

  // add the dimension for the cartesian_filter_
  StaticVector3 position_shape;
  if (cartesian_filter_) {
    position_shape = position_filter_->get_dimension();
    check_filter = false;
  } else {
    fatal_error("For the tally, " + name_ +
                ", the positional filter is required to do the legendre "
                "functional expansion tally.");
  }

  StaticVector6 tally_dimensions_(position_shape.begin(), position_shape.end());

  // add the dimension for energy_in_ only if greater than 1.
  if (energy_in_) {
    std::size_t ne = energy_in_->size();
    tally_dimensions_.insert(tally_dimensions_.begin(), ne);

    check_filter = check_filter && false;
  }

  if (check_filter == true) {
    fatal_error("for the " + name_ + " tally, no filter is provided.");
  }

  // Add the dimension for axis only when greater than 1
  std::size_t n_axis = axes.size();
  if (n_axis == 0) {
    fatal_error("For tally, " + name_ + " length of axis is zero.");
  }
  tally_dimensions_.push_back(n_axis);

  // Add the dimension for the order of FET only if greater than 0
  std::size_t fet_order_dimension = FET_order_ + 1;
  tally_dimensions_.push_back(fet_order_dimension);

  tally_avg.reallocate(tally_dimensions_);
  tally_avg.fill(0.0);

  tally_gen_score.reallocate(tally_dimensions_);
  tally_gen_score.fill(0.0);

  tally_var.reallocate(tally_dimensions_);
  tally_var.fill(0.0);
}

void LegendreFET::score_collision(const Particle& p, const Tracker& tktr,
                                  MaterialHelper& mat) {
  if (!(estimator_ == LegendreFET::Estimator::Collision)) {
    return;
  }

  StaticVector3 position_index = cartesian_filter_->get_indices(tktr);
  
  StaticVector6 indexes(position_index.begin(), position_index.end());

  std::size_t index_E;
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());
    
    if ( E_indx.has_value() )
      index_E = E_indx.value();

    else
      return;
  }

  // if the indexes_ is empty, then, we didn't get any scoring bin(s)
  if (indexes.empty()) {
    return;
  }

  const double Et = mat.Et(p.E());

  const double collision_score = particle_base_score(p, mat) / Et;
  
  // Check weather, given-axis size is more than 1.
  // if less than 1, there will not be any index
  const size_t axis_index = indexes.size();
  indexes.push_back(0);

  // Check weather, given- FET order is more than 1 or not.
  const size_t FET_index = indexes.size();
  indexes.push_back(0);

  // Variables for scoring
  double beta_n, scaled_loc;

  size_t it_axis = 0;
  for (auto& c : axes) {  // Loop over the different axis and indexing is done
                          // using the it_axis

    indexes[axis_index] = it_axis;

    switch (c) {
      case LegendreFET::Axis::X: {
        const double xmin_ = cartesian_filter_->x_min(position_index);
        const double xmax_ = cartesian_filter_->x_max(position_index);
        scaled_loc = 2. * (tktr.r().x() - xmin_) / (xmax_ - xmin_) - 1.;
      } break;

      case LegendreFET::Axis::Y: {
        const double ymin_ = cartesian_filter_->y_min(position_index);
        const double ymax_ = cartesian_filter_->y_max(position_index);
        scaled_loc = 2 * (tktr.r().y() - ymin_) / (ymax_ - ymin_) - 1;
      } break;

      case LegendreFET::Axis::Z: {
        const double zmin_ = cartesian_filter_->z_min(position_index);
        const double zmax_ = cartesian_filter_->z_max(position_index);
        scaled_loc = 2 * (tktr.r().z() - zmin_) / (zmax_ - zmin_) - 1;
      }
      default:
        break;
    }

    for (size_t i = 0; i < FET_order + 1;
         i++) {  // loop over differnt FET order
      beta_n =
          collision_score *
          legendre(i, scaled_loc);  // score for i-th order's basis function

      indexes[FET_index] = i;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score(indexes) += beta_n;
    }

    it_axis++;
  }
}

std::shared_ptr<LegendreFET> make_legendre_fet(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"]) {
    fatal_error("No valid name is provided.");
  }
  std::string legendre_fet_tally_name = node["name"].as<std::string>();

  // Check the estimator is given or not.
  if (!node["estimator"]) {
    fatal_error(
        "No, estimator is given for \"" + legendre_fet_tally_name +
        "\", therefore, collision estimator will be taken, by default.");
  }
  std::string estimator_name = node["estimator"].as<std::string>();

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"]) {
    fatal_error("No quantity is given for " + legendre_fet_tally_name +
                " tally.");
  }
  given_quantity = node["quantity"].as<std::string>();
  LegendreFET::Quantity quant;

  bool found_quantity = false;
  if (given_quantity == "flux") {
    quant = LegendreFET::Quantity::Flux;
    found_quantity = true;
  }

  if (given_quantity == "fission") {
    quant = LegendreFET::Quantity::Fission;
    found_quantity = true;
  }
  if (given_quantity == "absorption") {
    quant = LegendreFET::Quantity::Absorption;
    found_quantity = true;
  }
  if (given_quantity == "elastic") {
    quant = LegendreFET::Quantity::Elastic;
    found_quantity = true;
  }

  if (found_quantity == false) {
    fatal_error("For " + legendre_fet_tally_name +
                " tally, a unknown quantity is given.");
  }

  // Get the enrgy bounds, if any is given
  std::shared_ptr<EnergyFilter> energy_filter_ = nullptr;
  if (node["energy-bounds"]) {
    std::vector<double> energy_bounds =
        node["energy-bounds"].as<std::vector<double>>();
    energy_filter_ = std::make_shared<EnergyFilter>(energy_bounds);
  }

  // Get the cartesian type position filter
  std::shared_ptr<CartesianFilter> cartesian_filter_ =
      make_cartesian_filter(node);

  // Get the Legendre-FET order
  if (!node["FET-order"]) {
    fatal_error("Legendre-FET order is not given for the " +
                legendre_fet_tally_name + "tally.");
  }
  int FET_order_int_type = node["FET-order"].as<int>();
  if (FET_order_int_type < 0) {
    fatal_error("Legendre-FET order should not be negative.");
  }
  std::size_t FET_order_ = static_cast<std::size_t>(FET_order_int_type);

  // Get the Legendre-FET axes
  if (!node["FET-axis"]) {
    fatal_error("The axes for " + legendre_fet_tally_name +
                " should be given. ");
  }
  std::vector<std::string> fet_axis_type_string =
      node["FET-axis"].as<std::vector<std::string>>();

  std::vector<LegendreFET::Axis> fet_axis;
  fet_axis.reserve(fet_axis_type_string.size());

  for (auto& c : fet_axis_type_string) {
    if (c == "X") fet_axis.push_back(LegendreFET::Axis::X);

    if (c == "Y") fet_axis.push_back(LegendreFET::Axis::Y);

    if (c == "Z") fet_axis.push_back(LegendreFET::Axis::Z);
  }

  if (!(fet_axis.size() == fet_axis_type_string.size())) {
    fatal_error("A non-allowable name of the axis is given for " +
                legendre_fet_tally_name + ".");
  }

  // Make the Legendre FET tally class
  std::shared_ptr<LegendreFET> itally_legendre_fet_tally_ = nullptr;
  if (estimator_name == "collision") {
    itally_legendre_fet_tally_ = std::make_shared<LegendreFET>(
        cartesian_filter_, energy_filter_, fet_axis, FET_order_, quant,
        LegendreFET::Estimator::Collision, legendre_fet_tally_name);
  }

  if (estimator_name == "track-length") {
    /*itally_legendre_fet_tally_ = std::make_shared<LegendreFET>(FET_order_,
       LegendreFET::Quantity::Flux, LegendreFET::Estimator::TrackLength,
       legendre_fet_tally_name, fet_axis, cartesian_filter_, energy_filter_);
    */
    fatal_error(
        "The track-length for the legendre-FET is not avilable as asked in " +
        legendre_fet_tally_name + " tally. ");
  }

  if (itally_legendre_fet_tally_ == nullptr)
    fatal_error("Incorrect \"estimator\" is given.");

  return itally_legendre_fet_tally_;
}