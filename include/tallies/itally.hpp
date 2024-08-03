#ifndef ITALLY_H
#define ITALLY_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>

#include <yaml-cpp/yaml.h>
#include <xtensor/xarray.hpp>

#include <set>
#include <string>

struct Quantity {
  enum class Type {
    // Flux-like quantities
    Flux,
    RealFlux,
    ImagFlux,
    // Reaction rates
    Total,
    Fission,
    Absorption,
    Elastic,
    MT,
    // Heating
    Heating,
    // Source-like tallies
    Source,
    RealSource,
    ImagSource
  };

  Type type;
  std::uint32_t mt;

  bool operator==(const Quantity& R) const {
    if (this->type != R.type) return false;

    if (this->type == Type::MT) return this->mt == R.mt;

    return true;
  }

  bool operator!=(const Quantity& R) const { return !this->operator==(R); }
};

enum class Estimator { Collision, TrackLength, Source };

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

  // For source estimator
  virtual void score_source(const BankedParticle& p) = 0;

  // Record the avg and variance for the generation
  void record_generation(double mulitplier = 1.0);

  // evaluate tally
  virtual double evaluate(const Position& r, const double& E) const = 0;
  virtual std::vector<double> evaluate(
      const std::vector<std::pair<Position, double>> r_E) const = 0;

  // Clear the tally_gen_score
  void clear_generation() { tally_gen_score_.fill(0.0); }

  void set_net_weight(double weight_) {
    net_weight_ = weight_;
    inv_net_weight_ = 1.0 / weight_;
  }

  Estimator estimator() { return estimator_; }
  std::string estimator_str();

  const Quantity& quantity() { return quantity_; }
  std::string quantity_str();

  const std::string& name() const { return tally_name_; }

  virtual void write_tally() = 0;

  static const std::set<std::string> reserved_tally_names;

 protected:
  double particle_base_score(double E, double wgt, double wgt2,
                             MaterialHelper* mat) const;

  void var_to_std_on_mean();

  xt::xarray<double> tally_avg_;
  xt::xarray<double> tally_gen_score_;
  xt::xarray<double> tally_var_;

  size_t gen_ = 0;

  Quantity quantity_;
  Estimator estimator_;
  std::string tally_name_;

  double net_weight_, inv_net_weight_;
};

Quantity read_quantity(const YAML::Node& node, const std::string& name);

#endif