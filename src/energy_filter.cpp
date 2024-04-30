#include <tallies/energy_filter.hpp>

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node) {
  if (!node['energy-bounds']) {
    return nullptr;
  }
  std::vector<double> energy_bounds_ =
      node['energy-bounds'].as<std::vector<double>>();
  std::shared_ptr<EnergyFilter> energy_filter_ =
      std::make_shared<EnergyFilter>(energy_bounds_);

  return energy_filter_;
}
