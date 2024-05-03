#include <tallies/box_position_filter.hpp>
#include <tallies/general_tally.hpp>
#include <tallies/itally.hpp>
#include <tallies/legendre_fet.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <tallies/tallies.hpp>
#include <utils/output.hpp>

double ITally::particle_base_score(const Particle& p, MaterialHelper& mat) {
  double collision_score = inv_net_weight_;
  switch (quantity_) {
    case Quantity::Flux:
      collision_score *= p.wgt();
      break;

    case Quantity::Fission:
      collision_score *= p.wgt() * mat.Ef(p.E());
      break;

    case Quantity::Absorption:
      collision_score *= p.wgt() * mat.Ea(p.E());
      break;

    case Quantity::Elastic:
      collision_score *= p.wgt() * mat.Eelastic(p.E());
      break;
  }
  return collision_score;
}

void ITally::record_generation(double multiplier) {
  gen_++;

  const double dg = static_cast<double>(gen_);
  const double invs_dg = 1. / dg;
  const double invs_dg_m1 = 1. / (dg - 1.);  

  // All worker threads must send their generation score to the master.
  // Master must recieve all generations scores from workers and add
  // them to it's own generation score.
  mpi::Reduce_sum(tally_gen_score_.data_vector(), 0);

  // Only try to update average and variance is we are master, as worker
  // processes don't have copies of this data, so it will seg-fault.
  if (mpi::rank == 0) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < tally_avg_.size(); i++) {
      // Get new average

      double old_avg = tally_avg_[i];
      double val = tally_gen_score_[i] * multiplier;
      double avg = old_avg + (val - old_avg) * invs_dg;
      tally_avg_[i] = avg;

      // Get new variance
      if (gen_ > 1) {
        const double old_var = tally_var_[i];
        const double diff = val - old_avg;
        tally_var_[i] =
            old_var + (diff * diff * invs_dg) - (old_var * invs_dg_m1);
      }
    }
  }
  // Clear the entry for the tally_gen
  tally_gen_score_.fill(0.0);
}

std::string ITally::estimator_str() {
  switch (estimator_) {
    case Estimator::Collision:
      return "collision";
      break;

    case Estimator::TrackLength:
      return "track-length";
      break;

    default:
      return "unkown";
      break;
  }
}

std::string ITally::quantity_str() {
  switch (quantity_) {
    case Quantity::Flux:
      return "flux";
      break;
    case Quantity::Fission:
      return "fission";
      break;
    case Quantity::Absorption:
      return "absorption";
      break;
    case Quantity::Elastic:
      return "elastic";
      break;
    default:
      return "unknown";
      break;
  }
}

void make_itally(Tallies& tallies, const YAML::Node& node) {
  if (!node["name"]) {
    fatal_error("Tally name is not given.");
  }
  std::string tally_name_ = node["name"].as<std::string>();
  std::string tally_type_ = "general";
  if (!node["tally-type"]) {
    warning(
        "The tally-type is not given therefore, \"general\" will be assumed.");
  }
  tally_type_ = node["tally-type"].as<std::string>();

  std::shared_ptr<ITally> new_ITally = nullptr;
  if (tally_type_ == "general") {
    new_ITally = make_general_tally(node);
  }else if (tally_type_ == "legendre-fet") {
    new_ITally = make_legendre_fet(node);
  } else if (tally_type_ == "zernike-fet") {
    fatal_error("The zernike-fet for " + tally_name_ +
                " is not supported yet.");
  } else {
    fatal_error("For the " + tally_name_ +
                ", something went worng in the input.");
  }

  // Add the new_ITally of type ITally into the "tallies"
  tallies.add_ITally(new_ITally);
}