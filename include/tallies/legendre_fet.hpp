#ifndef LEGENDRE_FET_H
#define LEGENDRE_FET_H

#include <tallies/cartesian_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>
#include <utils/error.hpp>

#include <yaml-cpp/yaml.h>
#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<size_t, 6>;

#include <iostream>
#include <memory>
#include <vector>

class LegendreFET : public ITally {
 public:
  enum class Axis { X, Y, Z };

  LegendreFET(std::shared_ptr<CartesianFilter> position_filter,
              std::shared_ptr<EnergyFilter> energy_in,
              std::vector<LegendreFET::Axis> axes,
              std::vector<std::size_t> fet_order, Quantity quantity,
              Estimator estimator, std::string name);

  ~LegendreFET() = default;

  void score_collision(const Particle& p, const Tracker& trkr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& /*p*/, const Tracker& /*trkr*/,
                    double /*d_flight*/,
                    MaterialHelper& /*mat*/) override final {
    fatal_error("the track-length for the legendre-fet is not supoorted yet.");
  }

  void score_source(const BankedParticle& p) override final;

  double evaluate(const Position& r, const double& E) const override final;
  std::vector<double> evaluate(
      const std::vector<std::pair<Position, double>> r_E) const override final;

  void write_tally() override final;

  std::vector<std::size_t> get_fet_order() { return fet_order_; }

 private:
  std::shared_ptr<CartesianFilter> cartesian_filter_;
  std::shared_ptr<EnergyFilter> energy_in_;

  boost::container::static_vector<Axis, 3> axes_;

  // Note that the number of coefficient will be one more than the order.
  std::vector<std::size_t> fet_order_;
};

std::shared_ptr<LegendreFET> make_legendre_fet(const YAML::Node& node);

#endif