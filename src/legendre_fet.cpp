#include <tallies/itally.hpp>
#include <tallies/legendre_fet.hpp>
#include <tallies/tallies.hpp>
#include <utils/legendre.hpp>
#include <utils/output.hpp>

#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<size_t, 6>;

LegendreFET::LegendreFET(std::shared_ptr<CartesianFilter> position_filter,
                         std::shared_ptr<EnergyFilter> energy_in,
                         std::vector<LegendreFET::Axis> axes, size_t fet_order,
                         LegendreFET::Quantity quantity,
                         LegendreFET::Estimator estimator, std::string name)
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
    fatal_error("For the tally, " + name +
                ", the positional filter is required to do the legendre "
                "functional expansion tally.");
  }
  StaticVector3 position_shape = cartesian_filter_->get_shape();
  tally_shape.insert(tally_shape.end(), position_shape.begin(),
                     position_shape.end());

  // for legendre, check the size of axis vector should be between than 3.
  if (axes_.size() <= 1 && axes_.size() >= 3)
    fatal_error("the length of given axes for " + name +
                " tally is not between 1 to 3.");

  std::size_t n_axis = axes_.size();
  tally_shape.push_back(n_axis);

  // Add the dimension for the order of FET only if greater than 0
  std::size_t fet_order_dimension = fet_order_ + 1;
  tally_shape.push_back(fet_order_dimension);

  // reallocate and fill with zeros for the tally avg, gen-score and variance
  tally_avg_.reallocate(tally_shape);
  tally_avg_.fill(0.0);

  tally_gen_score_.reallocate(tally_shape);
  tally_gen_score_.fill(0.0);

  tally_var_.reallocate(tally_shape);
  tally_var_.fill(0.0);
}

void LegendreFET::score_collision(const Particle& p, const Tracker& tktr,
                                  MaterialHelper& mat) {
  StaticVector6 indices;
  // get the energy-index, if energy-filter exists
  if (energy_in_) {
    std::optional<std::size_t> E_indx = energy_in_->get_index(p.E());

    if (E_indx.has_value()) {
      std::size_t index_E = E_indx.value();
      indices.push_back(index_E);
    } else
      // Not inside any energy bin. Don't score.
      return;
  }

  // get the cartisian_filter indices
  StaticVector3 position_index = cartesian_filter_->get_indices(tktr);
  if (position_index.empty()) {
    // No bin is found, don't score.
    return;
  }

  indices.insert(indices.end(), position_index.begin(), position_index.end());

  const double Et = mat.Et(p.E());

  const double collision_score = particle_base_score(p, mat) / Et;

  // add one dimension for axis_index
  const size_t axis_index = indices.size();
  indices.push_back(0);

  // add one dimesnion for FET-index
  const size_t FET_index = indices.size();
  indices.push_back(0);

  // Variables for scoring
  double beta_n, scaled_loc;
  // Loop over the different axis and indexing is done
  for (std::size_t it_axis = 0; it_axis < axes_.size(); it_axis++) {
    // using the it_axis

    indices[axis_index] = it_axis;

    // get the scaled x, y, or z for legendre polynomial
    switch (axes_[it_axis]) {
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

    // loop over differnt FET order
    for (size_t i = 0; i < fet_order_ + 1; i++) {
      // score for i-th order's basis function
      beta_n = collision_score * legendre(i, scaled_loc);

      indices[FET_index] = i;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      tally_gen_score_(indices) += beta_n;
    }

    it_axis++;
  }
}

void LegendreFET::write_tally() {
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

std::shared_ptr<LegendreFET> make_legendre_fet(const YAML::Node& node) {
  // Check the name of the tally is given or not.
  if (!node["name"] && !node["name"].IsScalar()) {
    fatal_error("No valid name is provided.");
  }
  std::string legendre_fet_tally_name = node["name"].as<std::string>();

  // Check the estimator is given or not.
  if (!node["estimator"] && !node["estimator"].IsScalar()) {
    fatal_error("No valid estimator is given for " + legendre_fet_tally_name +
                " tally.");
  }
  std::string estimator_name = node["estimator"].as<std::string>();
  LegendreFET::Estimator estimator;
  if (estimator_name == "collision") {
    estimator = LegendreFET::Estimator::Collision;
  } else if (estimator_name == "track-length") {
    estimator = LegendreFET::Estimator::TrackLength;
    fatal_error(
        "The track-length for the legendre-fet is not supported as asked in " +
        legendre_fet_tally_name + " tally. ");
  } else {
    fatal_error("Incorrect \"estimator\" is given for the " +
                legendre_fet_tally_name + " tally.");
  }

  // Check for the quantity
  std::string given_quantity = "";
  if (!node["quantity"] && !node["quantity"].IsScalar()) {
    fatal_error("No valid quantity is given for " + legendre_fet_tally_name +
                " tally.");
  }
  given_quantity = node["quantity"].as<std::string>();
  LegendreFET::Quantity quant;

  if (given_quantity == "flux") {
    quant = LegendreFET::Quantity::Flux;
  } else if (given_quantity == "fission") {
    quant = LegendreFET::Quantity::Fission;
  } else if (given_quantity == "absorption") {
    quant = LegendreFET::Quantity::Absorption;
  } else if (given_quantity == "elastic") {
    quant = LegendreFET::Quantity::Elastic;
  } else {
    fatal_error("For " + legendre_fet_tally_name +
                " tally, a unknown quantity is given.");
  }

  // get the tallies instance
  auto& tallies = Tallies::instance();

  // Get the enrgy filter, if any is given
  std::shared_ptr<EnergyFilter> energy_filter = nullptr;
  if (node["energy-filters"]) {
    if (!node["energy-filters"].IsScalar()) {
      fatal_error("energy-filters is not provided as a scalar.");
    }
    std::size_t energy_id = node["energy-filters"].as<std::size_t>();
    energy_filter = tallies.get_energy_filter(energy_id);
    if (energy_filter == nullptr) {
      std::stringstream messg;
      messg << "for the tally " << legendre_fet_tally_name
            << ", the id: " << energy_id
            << ", for the energy-filters is not provided in the tally-filters.";
      fatal_error(messg.str());
    }
  }

  // Get the cartesian type position filter
  if (!node["position-filter-type"] ||
      !node["position-filter-type"].IsScalar()) {
    fatal_error("For " + legendre_fet_tally_name +
                ", a valid position-filter must be given.");
  }
  std::size_t position_id = node["position-filter-type"].as<std::size_t>();
  std::shared_ptr<CartesianFilter> cartesian_filter =
      tallies.get_cartesian_filter(position_id);
  if (cartesian_filter == nullptr) {
    std::stringstream messg;
    messg << "for tally " << legendre_fet_tally_name
          << ", the id: " << position_id
          << ", for the position-filter is not provided in the tally-filters.";
    fatal_error(messg.str());
  }

  // Get the legendre-fet order
  if (!node["order"] && !node["order"].IsScalar()) {
    fatal_error("legendre-fet order is not given for the " +
                legendre_fet_tally_name + "tally.");
  }
  int fet_order_int_type = node["order"].as<int>();
  if (fet_order_int_type < 0) {
    fatal_error("legendre-fet order should not be negative.");
  }
  std::size_t FET_order = static_cast<std::size_t>(fet_order_int_type);

  // Get the legendre-fet axes
  if (!node["axes"] && !node["axes"].IsSequence()) {
    fatal_error("A valid \"axes\" list for " + legendre_fet_tally_name +
                " should be given. ");
  }
  std::vector<std::string> fet_axis_type_string =
      node["axes"].as<std::vector<std::string>>();

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
  std::shared_ptr<LegendreFET> itally_legendre_fet_tally_ =
      std::make_shared<LegendreFET>(cartesian_filter, energy_filter, fet_axis,
                                    FET_order, quant, estimator,
                                    legendre_fet_tally_name);

  return itally_legendre_fet_tally_;
}