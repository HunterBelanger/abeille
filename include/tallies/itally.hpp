#ifndef ITALLY_H
#define ITALLY_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>

#include <ndarray.hpp>

#include <set>
#include <string>

enum class Quantity {
  Flux,
  Total,
  Fission,
  Absorption,
  Elastic,
  RealFlux,
  ImagFlux
};

enum class Estimator { Collision, TrackLength };

class ITally {
 public:
  ITally(Quantity quantity, Estimator estimator, std::string name);
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

  const std::string& name() const { return tally_name_; }

  virtual void write_tally() = 0;

  static const std::set<std::string> reserved_tally_names;

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

Quantity read_quantity(const std::string& quant_str,
                       const std::string& tally_name);

#endif