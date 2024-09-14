#include <tallies/energy_filter.hpp>
#include <utils/error.hpp>

#include <sstream>

EnergyFilter::EnergyFilter(const std::vector<double> energy_bounds,
                           std::size_t id)
    : energy_bounds_(energy_bounds), id_(id) {
  // Check the size of energy_bounds should be multiple of 2
  // since any bounds has 2 end-point-boundary-elements to define the range.
  if (energy_bounds_.size() < 2) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id_
         << " has less than 2 energy bounds.";
    fatal_error(mssg.str());
  }

  if (std::is_sorted(energy_bounds_.begin(), energy_bounds_.end()) == false) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id_ << " has unsorted energy bounds.";
    fatal_error(mssg.str());
  }

  if (energy_bounds_[0] < 0.0) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id_ << " has negative energy bounds.";
    fatal_error(mssg.str());
  }
}

std::optional<std::size_t> EnergyFilter::get_index(const double& E) const {
  for (std::size_t i = 0; i < energy_bounds_.size() - 1; i++) {
    if (energy_bounds_[i] <= E && E <= energy_bounds_[i + 1]) {
      return i;
    }
  }
  // if we get here, the energy index is not found.
  return std::nullopt;
}

void EnergyFilter::write_to_hdf5(H5::Group& grp) const {
  // Save id in attributes
  if (grp.hasAttribute("id")) {
    grp.deleteAttribute("id");
  }
  grp.createAttribute("id", this->id());

  // Save energy bounds in a dataset
  if (grp.exist("energy-bounds")) {
    grp.unlink("energy-bounds");
  }
  grp.createDataSet("energy-bounds", this->energy_bounds());
}

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node) {
  if (!node["id"] || !node["id"].IsScalar()) {
    fatal_error("invalid id is given for energy-fitler.");
  }
  std::size_t id = node["id"].as<std::size_t>();

  if (!node["energy-bounds"].IsSequence() ||
      (node["energy-bounds"].size() < 2)) {
    std::stringstream mssg;
    mssg << "Invalid enery-bounds are given on energy-filter with id " << id
         << ".";
    fatal_error(mssg.str());
  }

  std::vector<double> energy_bounds =
      node["energy-bounds"].as<std::vector<double>>();

  if (energy_bounds.size() < 2) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id << " has less than 2 energy bounds.";
    fatal_error(mssg.str());
  }

  if (std::is_sorted(energy_bounds.begin(), energy_bounds.end()) == false) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id << " has unsorted energy bounds.";
    fatal_error(mssg.str());
  }

  if (energy_bounds[0] < 0.0) {
    std::stringstream mssg;
    mssg << "Energy filter with id " << id << " has negative energy bounds.";
    fatal_error(mssg.str());
  }

  return std::make_shared<EnergyFilter>(energy_bounds, id);
}
