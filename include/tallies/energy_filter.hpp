#ifndef ENREGY_FILTER_H
#define ENREGY_FILTER_H

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <string>
#include <vector>
#include <optional>


class EnergyFilter {
 public:
  EnergyFilter(const std::vector<double> energy_bounds);

  ~EnergyFilter() = default;

  std::optional<std::size_t> get_index(const double& E) const;

  std::size_t size() { return energy_bounds_.size() - 1; }
  const std::vector<double>& energy_bounds() const { return energy_bounds_; }

  std::string type_str() { return "energyfilter"; }

 private:
  std::vector<double> energy_bounds_;
};

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node);

#endif