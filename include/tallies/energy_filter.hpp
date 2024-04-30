#ifndef ENREGY_FILTER_H
#define ENREGY_FILTER_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include <tallies/position_filter.hpp>
#include <utils/error.hpp>

class EnergyFilter {
 public:
  EnergyFilter(const std::vector<double> energy_bounds)
      : energy_bounds_(energy_bounds) {
    // Check the size of energy_bounds should be multiple of 2
    // since any bounds has 2 end-point-boundary-elements to define the range.
    if (energy_bounds.size() < 2) {
      fatal_error("Size of energy bounds is less than 2 in EnergyFilter.");
    }

    if (std::is_sorted(energy_bounds_.begin(), energy_bounds_.end()) == false) {
      fatal_error("Energy bounds must be sorted on EnergyFilter.");
    }

    if (energy_bounds_[0] < 0.0) {
      fatal_error(
          "First element of energy bounds is less than zero in EnergyFilter.");
    }
  }

  ~EnergyFilter() = default;

  bool get_index(const double& E, std::size_t& index_E) const {
    // index_E = -1;
    for (std::size_t i = 0; i < energy_bounds_.size() - 1; i++) {
      if (E >= energy_bounds_[i] && E <= energy_bounds_[i + 1]) {
        index_E = i;
        return true;
        break;
      }
    }
    return false;
  }

  std::size_t size() { return energy_bounds_.size() - 1; }
  const std::vector<double>& energy_bounds() const { return energy_bounds_; }

  // Perhaps Not Needed
  FilterType type() { return FilterType::Energy_Filter; }
  std::string type_str() { return "energyfilter"; }

 private:
  std::vector<double> energy_bounds_;
};

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node);

#endif