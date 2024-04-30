#ifndef ITALLY_H
#define ITALLY_H

#include <memory>

#include <yaml-cpp/yaml.h>
#include <ndarray.hpp>

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <tallies/box_position_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/mesh_position_filter.hpp>
#include <tallies/position_filter.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/position.hpp>

#include <boost/container/static_vector.hpp>

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

  ITally(Quantity quantity, Estimator estimator, std::string name_)
      : quantity_(quantity), estimator_(estimator), tally_name(name_) {}

  virtual ~ITally() = default;

  // For collision estmator
  virtual void score_collision(const Particle& p, const Tracker& tktr,
                               MaterialHelper& mat) = 0;

  // For track-length estimator
  // virtual void score_flight(const Particle& p, double d_flight
  // ,MaterialHelper& mat ) = 0;
  virtual void score_flight(const Particle& p, const Tracker& trkr,
                            double d_flight, MaterialHelper& mat) = 0;

  // Record the avg and variance for the generation
  void record_generation(double mulitplier = 1.0);

  void set_net_weight(double weight_) { net_weight_ = weight_; }

  std::string estimator_str();
  std::string quantity_str();

  void write_tally();

 protected:
  double particle_base_score(const Particle& p, MaterialHelper& mat);

 protected:
  NDArray<double> tally_avg;
  NDArray<double> tally_gen_score;
  NDArray<double> tally_var;

  size_t gen_ = 0;

  Quantity quantity_;
  Estimator estimator_;
  std::string tally_name;

  double net_weight_;
};

#endif