#include <tallies/general_tally.hpp>
#include <tallies/itally.hpp>
#include <tallies/legendre_fet.hpp>
#include <tallies/regular_cartesian_mesh_filter.hpp>
#include <tallies/tallies.hpp>
#include <tallies/zernike_fet.hpp>
#include <utils/mpi.hpp>

const std::set<std::string> ITally::reserved_tally_names{
    "families",
    "pair-dist-sqrd",
    "entropy",
    "total-pre-cancel-entropy",
    "neg-pre-cancel-entropy",
    "pos-pre-cancel-entropy",
    "total-post-cancel-entropy",
    "neg-post-cancel-entropy",
    "pos-post-cancel-entropy",
    "empty-entropy-frac",
    "Nnet",
    "Ntot",
    "Npos",
    "Nneg",
    "Wnet",
    "Wtot",
    "Wpos",
    "Wneg",
    "kcol",
    "ktrk",
    "kabs",
    "leakage",
    "mig-area"};

ITally::ITally(Quantity quantity, Estimator estimator, std::string name)
    : quantity_(quantity), estimator_(estimator), tally_name_(name) {
  if (reserved_tally_names.contains(name)) {
    fatal_error("The tally name " + name + " is reserved.");
  }
}

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

    case Quantity::Total:
      collision_score *= p.wgt() * mat.Et(p.E());
      break;

    case Quantity::RealFlux:
      collision_score *= p.wgt();
      break;  

    case Quantity::ImgFlux:
      collision_score *= p.wgt2();
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
    for (std::size_t i = 0; i < tally_avg_.size(); i++) {
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
  }

  // NEVER GETS HERE
  return "unknown";
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
    case Quantity::Total:
      return "total";
      break;
    case Quantity::RealFlux:
      return "real-flux";
      break;
    case Quantity::ImgFlux:
      return "imarginary-flux";
      break;
  }

  // NEVER GETS HERE
  return "unknown";
}

void make_itally(Tallies& tallies, const YAML::Node& node) {
  if (!node["name"] || node["name"].IsScalar() == false) {
    fatal_error("Tally name is not given.");
  }
  std::string tally_name = node["name"].as<std::string>();
  if (ITally::reserved_tally_names.contains(tally_name)) {
    fatal_error("The tally name " + tally_name + " is reserved.");
  }

  std::string tally_type = "general";
  if (node["tally-type"] && node["tally-type"].IsScalar()) {
    tally_type = node["tally-type"].as<std::string>();
  } else if (node["tally-type"]) {
    fatal_error("Tally " + tally_name + " had invalid type entry.");
  }

  std::shared_ptr<ITally> t = nullptr;
  if (tally_type == "general") {
    t = make_general_tally(node);
  } else if (tally_type == "legendre-fet") {
    t = make_legendre_fet(node);
  } else if (tally_type == "zernike-fet") {
    t = make_zernike_fet(node);
  } else {
    fatal_error("Unknown tally type " + tally_type + " found in tally " +
                tally_name + ".");
  }

  // Add the new_ITally of type ITally into the "tallies"
  tallies.add_ITally(t);
}