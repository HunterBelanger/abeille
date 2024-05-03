#include <tallies/energy_filter.hpp>
#include <utils/error.hpp>

EnergyFilter::EnergyFilter(const std::vector<double> energy_bounds)
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

std::optional<std::size_t> EnergyFilter::get_index(const double& E) const {
  for (std::size_t i = 0; i < energy_bounds_.size() - 1; i++) {
    if (energy_bounds_[i] <= E && E <= energy_bounds_[i + 1]) {
      return i;
      break;
    }
  }
  // if we get here, the energy index is not found.
  return std::nullopt;
}

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node) {
  if (!node["energy-bounds"].IsSequence() ||
      (node["energy-bounds"].size() < 2)) {
    fatal_error("Invalid enery-bounds are given.");
  }
  std::vector<double> energy_bounds_ =
      node["energy-bounds"].as<std::vector<double>>();

  return std::make_shared<EnergyFilter>(energy_bounds_);
}
