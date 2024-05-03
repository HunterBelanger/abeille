#ifndef GENERAL_TALLY_H
#define GENERAL_TALLY_H

#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>
#include <tallies/position_filter.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class GeneralTally : public ITally {
 public:
  GeneralTally(std::shared_ptr<PositionFilter> position_filter,
               std::shared_ptr<EnergyFilter> energy_in, Quantity quantity,
               Estimator estimator, std::string name_);

  ~GeneralTally() = default;

  void score_collision(const Particle& p, const Tracker& tktr,
                       MaterialHelper& mat) override final;

  void score_flight(const Particle& p, const Tracker& trkr, double d_flight,
                    MaterialHelper& mat) override final;

 private:
  std::shared_ptr<PositionFilter> position_filter_;
  std::shared_ptr<EnergyFilter> energy_in_;
};

std::shared_ptr<GeneralTally> make_general_tally(const YAML::Node& node);

#endif