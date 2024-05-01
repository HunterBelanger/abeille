#ifndef LEGENDRE_FET_H
#define LEGENDRE_FET_H

#include <tallies/itally.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>

#include <yaml-cpp/yaml.h>
#include <ndarray.hpp>
#include <boost/container/static_vector.hpp>
using StaticVector6 = boost::container::static_vector<size_t, 6>;

#include <iostream>
#include <memory>
#include <vector>

class LegendreFET : public ITally {
 public:
  enum class Axis { X, Y, Z };

  LegendreFET(std::shared_ptr<CartesianFilter> position_filter_,
              std::shared_ptr<EnergyFilter> energy_in,
              std::vector<LegendreFET::Axis> axes_, size_t FET_order_,
              LegendreFET::Quantity quantity_,
              LegendreFET::Estimator estimator_, std::string name_);

  ~LegendreFET() = default;

  void score_collision(const Particle& p, const Tracker& tktr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& /*p*/, const Tracker& /*trkr*/,
                    double /*d_flight*/,
                    MaterialHelper& /*mat*/) override final {}

  size_t get_fet_order() { return FET_order; }

 private:
  std::shared_ptr<CartesianFilter> cartesian_filter_;
  std::shared_ptr<EnergyFilter> energy_in_;
  std::shared_ptr<EnergyFilter> energy_out_ = nullptr;

  boost::container::static_vector<Axis, 3> axes;
  size_t FET_order;  // Note that the number of coefficient will be one more
                     // than the order.
};

std::shared_ptr<LegendreFET> make_legendre_fet(const YAML::Node& node);

#endif