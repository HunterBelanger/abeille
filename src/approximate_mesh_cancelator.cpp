/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#include <simulation/approximate_mesh_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <ndarray.hpp>

#include <algorithm>
#include <cmath>

ApproximateMeshCancelator::ApproximateMeshCancelator(Position low, Position hi,
                                                     uint32_t Nx, uint32_t Ny,
                                                     uint32_t Nz)
    : r_low(low),
      r_hi(hi),
      energy_edges(),
      shape{Nx, Ny, Nz, 1},
      dx(0.),
      dy(0.),
      dz(0.),
      bins() {
  // Make sure the low points are all lower than the high points
  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    fatal_error(
        "Low position is not lower than hi position in "
        "ApproximateMeshCancelator.");
  }

  dx = (r_hi.x() - r_low.x()) / static_cast<double>(shape[0]);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(shape[1]);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(shape[2]);
}

ApproximateMeshCancelator::ApproximateMeshCancelator(
    Position low, Position hi, uint32_t Nx, uint32_t Ny, uint32_t Nz,
    std::vector<double> energy_bounds)
    : r_low(low),
      r_hi(hi),
      energy_edges(energy_bounds),
      shape{Nx, Ny, Nz, 1},
      dx(0.),
      dy(0.),
      dz(0.),
      bins() {
  // Make sure the low points are all lower than the high points
  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    fatal_error(
        "Low position is not lower than hi position in "
        "ApproximateMeshCancelator.");
  }

  dx = (r_hi.x() - r_low.x()) / static_cast<double>(shape[0]);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(shape[1]);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(shape[2]);

  // Make sure energy bins are valid
  if (energy_edges.size() < 2) {
    fatal_error(
        "energy_edges must have at least two entries in "
        "ApproximateMeshCancelator.");
  }

  if (!std::is_sorted(energy_edges.begin(), energy_edges.end())) {
    fatal_error("energy_edges must be sorted in ApproximateMeshCancelator.");
  }

  if (energy_edges.front() < 0.) {
    fatal_error(
        "All energy_edges must be greater than or equal to zero in "
        "ApproximateMeshCancelator.");
  }

  shape[3] = static_cast<uint32_t>(energy_edges.size() - 1);
}

bool ApproximateMeshCancelator::add_particle(BankedParticle& p) {
  // Get bin indicies for spacial coordinates
  int i = static_cast<int>(std::floor((p.r.x() - r_low.x()) / dx));
  int j = static_cast<int>(std::floor((p.r.y() - r_low.y()) / dy));
  int k = static_cast<int>(std::floor((p.r.z() - r_low.z()) / dz));

  // Get energy index with linear search
  int l = -1;
  if (energy_edges.empty() == false) {
    for (size_t e = 0; e < energy_edges.size() - 1; e++) {
      if (energy_edges[e] <= p.E && p.E <= energy_edges[e + 1]) {
        l = static_cast<int>(e);
        break;
      }
    }
  } else {
    l = 0;
  }

  // If one of the indices is less than zero, don't keep the particle
  if (i < 0 || j < 0 || k < 0 || l < 0) {
    return false;
  }

  // if one of the indices is too large, don't keep the particle
  if (i >= static_cast<int>(shape[0]) || j >= static_cast<int>(shape[1]) ||
      k >= static_cast<int>(shape[2]) || l >= static_cast<int>(shape[3])) {
    return false;
  }

  // The particle fits into the mesh somewhere
  auto key = [shape = shape](const int& i, const int& j, const int& k,
                             const int& l) {
    return l + static_cast<int>(shape[3]) *
                   (k + static_cast<int>(shape[2]) *
                            (j + static_cast<int>(shape[1]) * i));
  };

  int bin_key = key(i, j, k, l);

  if (bins.find(bin_key) == bins.end()) {
    bins[bin_key] = std::vector<BankedParticle*>();
  }

  bins[bin_key].push_back(&p);

  return true;
}

std::vector<int> ApproximateMeshCancelator::sync_keys() {
  // Each node collects all keys
  std::vector<int> keys;
  keys.reserve(bins.size());
  for (auto& key_bin_pair : bins) {
    keys.push_back(key_bin_pair.first);
  }
  std::set<int> key_set;

  // Put master keys into the keyset
  if (mpi::rank == 0) {
    std::copy(keys.begin(), keys.end(), std::inserter(key_set, key_set.end()));
    keys.clear();
  }

  // For every Node starting at 1, send its keys to master and add to key_set
  for (int i = 1; i < mpi::size; i++) {
    if (mpi::rank == i) {
      mpi::Send(keys, 0);
      keys.clear();
    } else if (mpi::rank == 0) {
      mpi::Recv(keys, i);
      std::copy(keys.begin(), keys.end(),
                std::inserter(key_set, key_set.end()));
    }
  }

  // Move keys from key_set back to keys in master
  if (mpi::rank == 0) {
    keys.clear();
    keys.assign(key_set.begin(), key_set.end());
    key_set.clear();
  }

  // Send keys to all nodes from Master
  mpi::Bcast(keys, 0);

  return keys;
}

void ApproximateMeshCancelator::perform_cancellation_loop(pcg32& /*rng*/) {
  // Get keys of all non empty bins
  std::vector<int> keys = sync_keys();

  for (const auto key : keys) {
    std::uint64_t n_total = 0;
    double sum_wgt = 0.;
    double sum_wgt2 = 0.;

    // Go through all particles in the bin and add up positives and negatives
    // and get total sum
    auto& bin = bins[key];
    for (const auto& p : bin) {
      n_total++;
      sum_wgt += p->wgt;
      sum_wgt2 += p->wgt2;
    }
    // Sum weights of all particles in each bin across all nodes
    mpi::Allreduce_sum(sum_wgt);
    mpi::Allreduce_sum(sum_wgt2);

    // Sum total number of positive and negative particles across all nodes
    mpi::Allreduce_sum(n_total);

    // Set the avg weights
    const double avg_wgt = sum_wgt / n_total;
    const double avg_wgt2 = sum_wgt2 / n_total;

    // Loop through particles in the bin once again and perform cancellation if
    // necessary
    for (const auto& p : bin) {
       p->wgt = avg_wgt;
       p->wgt2 = avg_wgt2;
    }
    bins[key].clear();
  }
  keys.clear();
}
void ApproximateMeshCancelator::perform_cancellation_vector(pcg32& /*rng*/) {
  // Get keys of all non empty bins
  std::vector<int> keys = sync_keys();

  // change to uint16
  NDArray<double> wgts({2, keys.size()});
  std::vector<uint16_t> n_totals(keys.size(), 0);
  
  for (std::size_t i = 0; i < keys.size(); i++) {
    const auto key = keys[i];
    std::uint64_t n_total = 0;
    double sum_wgt = 0.;
    double sum_wgt2 = 0.;

    // Go through all particles in the bin and add up positives and negatives
    // and get total sum
    for (const auto& p : bins[key]) {
      n_total++;
      sum_wgt += p->wgt;
      sum_wgt2 += p->wgt2;
    }

    // Push the counts to the vectors
    n_totals[i] = n_total;
    wgts(0, i) = sum_wgt;
    wgts(1, i) = sum_wgt2;
  }

  // Sum the vectors across all nodes
  mpi::Allreduce_sum(wgts.data_vector());
  mpi::Allreduce_sum(n_totals);

  // all the vectors have size keys.size() so we use variable x to index them
  // since they should match to keys
  for (std::size_t i = 0; i < keys.size();i++) {
    const auto key = keys[i];

    // Set the avg weights
    const double n = static_cast<double>(n_totals[i]);
    const double avg_wgt = wgts(0,i) / n;
    const double avg_wgt2 = wgts(1,i) / n;

    for (auto& p : bins[key]) {
       p->wgt = avg_wgt;
       p->wgt2 = avg_wgt2;
    }

    bins[key].clear();
  }

  keys.clear();
}

void ApproximateMeshCancelator::perform_cancellation(pcg32& /*rng*/) {
  // Go through all bins in the mesh
  for (auto& key_bin_pair : bins) {
    auto& bin = key_bin_pair.second;

    // Only do cancelation if we have more than one particle per bin
    if (bin.size() > 1) {
      // Only do cancellation if we have differing signs
      bool has_pos_w1 = false;
      bool has_neg_w1 = false;
      bool has_pos_w2 = false;
      bool has_neg_w2 = false;
      double sum_wgt = 0.;
      double sum_wgt2 = 0.;

      // Go through all particles in the bin
      for (const auto& p : bin) {
        if (p->wgt > 0.)
          has_pos_w1 = true;
        else if (p->wgt < 0.)
          has_neg_w1 = true;
        sum_wgt += p->wgt;

        if (p->wgt2 > 0.)
          has_pos_w2 = true;
        else if (p->wgt2 < 0.)
          has_neg_w2 = true;
        sum_wgt2 += p->wgt2;
      }

      // Get average weights
      double N = static_cast<double>(bin.size());
      double avg_wgt = sum_wgt / N;
      double avg_wgt2 = sum_wgt2 / N;

      // Go through all particles and change their weights
      for (auto& p : bin) {
        if (has_pos_w1 && has_neg_w1) p->wgt = avg_wgt;
        if (has_pos_w2 && has_neg_w2) p->wgt2 = avg_wgt2;
      }
    }

    bin.clear();
  }
}

std::vector<BankedParticle> ApproximateMeshCancelator::get_new_particles(
    pcg32& /*rng*/) {
  return {};
}

void ApproximateMeshCancelator::clear() { bins.clear(); }

std::shared_ptr<ApproximateMeshCancelator> make_approximate_mesh_cancelator(
    const YAML::Node& node) {
  // Get low
  if (!node["low"] || !node["low"].IsSequence() || !(node["low"].size() == 3)) {
    fatal_error("No valid low entry for approximate mesh cancelator.");
  }

  double xl = node["low"][0].as<double>();
  double yl = node["low"][1].as<double>();
  double zl = node["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!node["hi"] || !node["hi"].IsSequence() || !(node["hi"].size() == 3)) {
    fatal_error("No valid hi entry for approximate mesh cancelator.");
  }

  double xh = node["hi"][0].as<double>();
  double yh = node["hi"][1].as<double>();
  double zh = node["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get shape
  if (!node["shape"] || !node["shape"].IsSequence() ||
      !(node["shape"].size() == 3)) {
    fatal_error("No valid shape entry for approximate mesh cancelator.");
  }

  uint32_t Nx = node["shape"][0].as<uint32_t>();
  uint32_t Ny = node["shape"][1].as<uint32_t>();
  uint32_t Nz = node["shape"][2].as<uint32_t>();

  if (!node["energy-bounds"]) {
    return std::make_shared<ApproximateMeshCancelator>(r_low, r_hi, Nx, Ny, Nz);
  }

  if (!node["energy-bounds"].IsSequence()) {
    fatal_error(
        "No valid energy-bounds entry for approximate mesh cancelator.");
  }

  std::vector<double> energy_bounds =
      node["energy-bounds"].as<std::vector<double>>();

  Output::instance().write(" Using ApproximateMeshCancelator.\n");

  return std::make_shared<ApproximateMeshCancelator>(r_low, r_hi, Nx, Ny, Nz,
                                                     energy_bounds);
}
