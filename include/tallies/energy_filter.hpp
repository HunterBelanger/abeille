#ifndef ENREGY_FILTER_H
#define ENREGY_FILTER_H

#include <yaml-cpp/yaml.h>

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <iostream>
#include <optional>
#include <string>
#include <vector>

class EnergyFilter {
 public:
  EnergyFilter(const std::vector<double> energy_bounds, std::size_t id);

  ~EnergyFilter() = default;

  std::optional<std::size_t> get_index(const double& E) const;

  std::size_t size() const { return energy_bounds_.size() - 1; }
  const std::vector<double>& energy_bounds() const { return energy_bounds_; }

  std::string type_str() const { return "energy-filter"; }

  std::size_t id() const { return id_; }

  void write_to_hdf5(H5::Group& grp) const;

 private:
  std::vector<double> energy_bounds_;
  std::size_t id_;
};

std::shared_ptr<EnergyFilter> make_energy_filter(const YAML::Node& node);

#endif