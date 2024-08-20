#include <tallies/itally.hpp>
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

double ITally::particle_base_score(double E, double wgt, double wgt2,
                                   MaterialHelper* mat) const {
  double collision_score = inv_net_weight_;
  switch (quantity_.type) {
    case Quantity::Type::Flux:
      collision_score *= wgt;
      break;

    case Quantity::Type::RealFlux:
      collision_score *= wgt;
      break;

    case Quantity::Type::ImagFlux:
      collision_score *= wgt2;
      break;

    case Quantity::Type::Total:
      collision_score *= wgt * mat->Et(E);
      break;

    case Quantity::Type::Fission:
      collision_score *= wgt * mat->Ef(E);
      break;

    case Quantity::Type::Absorption:
      collision_score *= wgt * mat->Ea(E);
      break;

    case Quantity::Type::Elastic:
      collision_score *= wgt * mat->Eelastic(E);
      break;

    case Quantity::Type::MT: {
      const double Emt = mat->Emt(quantity_.mt, E);
      collision_score *= wgt * Emt;
    } break;

    case Quantity::Type::Heating: {
      const double heating = mat->heating(E);
      collision_score *= wgt * heating;
    } break;

    case Quantity::Type::Source:
      collision_score *= wgt;
      break;

    case Quantity::Type::RealSource:
      collision_score *= wgt;
      break;

    case Quantity::Type::ImagSource:
      collision_score *= wgt2;
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
  std::span<double> gen_vals(tally_gen_score_.data(), tally_gen_score_.size());
  mpi::Reduce_sum(gen_vals, 0);

  // Only try to update average and variance is we are master, as worker
  // processes don't have copies of this data, so it will seg-fault.
  if (mpi::rank == 0) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < tally_avg_.size(); i++) {
      // Get new average

      double old_avg = tally_avg_.flat(i);
      double val = tally_gen_score_.flat(i) * multiplier;
      double avg = old_avg + (val - old_avg) * invs_dg;
      tally_avg_.flat(i) = avg;

      // Get new variance
      if (gen_ > 1) {
        const double old_var = tally_var_.flat(i);
        const double diff = val - old_avg;
        tally_var_.flat(i) =
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

    case Estimator::Source:
      return "source";
      break;
  }

  // NEVER GETS HERE
  return "unknown";
}

std::string ITally::quantity_str() {
  switch (quantity_.type) {
    case Quantity::Type::Flux:
      return "flux";
      break;
    case Quantity::Type::RealFlux:
      return "real-flux";
      break;
    case Quantity::Type::ImagFlux:
      return "imag-flux";
      break;
    case Quantity::Type::Total:
      return "total";
      break;
    case Quantity::Type::Fission:
      return "fission";
      break;
    case Quantity::Type::Absorption:
      return "absorption";
      break;
    case Quantity::Type::Elastic:
      return "elastic";
      break;
    case Quantity::Type::MT:
      return "mt";
      break;
    case Quantity::Type::Heating:
      return "heating";
      break;
    case Quantity::Type::Source:
      return "source";
      break;
    case Quantity::Type::RealSource:
      return "real-source";
      break;
    case Quantity::Type::ImagSource:
      return "imag-source";
      break;
  }

  // NEVER GETS HERE
  return "unknown";
}

void ITally::var_to_std_on_mean() {
  const double invs_g = 1. / static_cast<double>(gen_);
  tally_var_ *= invs_g;
  tally_var_ = xt::sqrt(tally_var_);
}

Quantity read_quantity(const YAML::Node& node, const std::string& name) {
  if (!node["quantity"] || node["quantity"].IsScalar() == false) {
    std::stringstream mssg;
    mssg << "Tally " << name << " is missing a valid quantity entry.";
    fatal_error(mssg.str());
  }
  std::string quant_str = node["quantity"].as<std::string>();

  if (quant_str == "flux") {
    return {Quantity::Type::Flux, 0};
  } else if (quant_str == "fission") {
    return {Quantity::Type::Fission, 0};
  } else if (quant_str == "absorption") {
    return {Quantity::Type::Absorption, 0};
  } else if (quant_str == "elastic") {
    return {Quantity::Type::Elastic, 0};
  } else if (quant_str == "total") {
    return {Quantity::Type::Total, 0};
  } else if (quant_str == "real-flux") {
    return {Quantity::Type::RealFlux, 0};
  } else if (quant_str == "imag-flux") {
    return {Quantity::Type::ImagFlux, 0};
  } else if (quant_str == "mt") {
    Quantity quant{Quantity::Type::MT, 0};

    if (settings::energy_mode == settings::EnergyMode::MG) {
      // Can't do an MT tally in MG mode !
      fatal_error("Cannot do MT tallies in multi-group mode.");
    }

    // Check for mt
    if (!node["mt"] || node["mt"].IsScalar() == false) {
      std::stringstream mssg;
      mssg << "Tally " << name
           << " has a quantity of mt, but not provided mt value.";
      fatal_error(mssg.str());
    }

    std::uint32_t tmp_mt = node["mt"].as<std::uint32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      std::stringstream mssg;
      mssg << "Tally " << name << " has an invalid mt value " << tmp_mt << ".";
      fatal_error(mssg.str());
    }

    quant.mt = tmp_mt;

    return quant;
  } else if (quant_str == "heating") {
    return {Quantity::Type::Heating, 0};
  } else if (quant_str == "source") {
    return {Quantity::Type::Source, 0};
  } else if (quant_str == "real-source") {
    return {Quantity::Type::RealSource, 0};
  } else if (quant_str == "imag-source") {
    return {Quantity::Type::ImagSource, 0};
  } else {
    fatal_error("Tally " + name + " has unknown quantity \"" + quant_str +
                "\".");
  }

  // SHOULD NEVER GET HERE
  return {Quantity::Type::Flux, 0};
}