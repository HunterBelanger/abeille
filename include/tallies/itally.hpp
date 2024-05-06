#ifndef ITALLY_H
#define ITALLY_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>

#include <ndarray.hpp>

const static std::set<std::string> reserved_tally_names{
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

class ITally {
 public:
  enum class Quantity {
    Flux,
    Fission,
    Absorption,
    Elastic,
  };

  enum class Estimator { Collision, TrackLength };

  ITally() = default;

  ITally(Quantity quantity, Estimator estimator, std::string name)
      : quantity_(quantity), estimator_(estimator), tally_name_(name) {
    if (reserved_tally_names.contains(name)) {
      fatal_error("The tally name " + name + " is reserved.");
    }
  }

  virtual ~ITally() = default;

  // For collision estmator
  virtual void score_collision(const Particle& p, const Tracker& tktr,
                               MaterialHelper& mat) = 0;

  // For track-length estimator
  virtual void score_flight(const Particle& p, const Tracker& trkr,
                            double d_flight, MaterialHelper& mat) = 0;

  // Record the avg and variance for the generation
  void record_generation(double mulitplier = 1.0);

  // Clear the tally_gen_score
  void clear_generation() { tally_gen_score_.fill(0.0); }

  void set_net_weight(double weight_) {
    net_weight_ = weight_;
    inv_net_weight_ = 1.0 / weight_;
  }

  Estimator estimator() { return estimator_; }
  std::string estimator_str();

  Quantity quantity() { return quantity_; }
  std::string quantity_str();

  virtual void write_tally() = 0;

 protected:
  double particle_base_score(const Particle& p, MaterialHelper& mat);

  NDArray<double> tally_avg_;
  NDArray<double> tally_gen_score_;
  NDArray<double> tally_var_;

  size_t gen_ = 0;

  Quantity quantity_;
  Estimator estimator_;
  std::string tally_name_;

  double net_weight_, inv_net_weight_;
};

#endif