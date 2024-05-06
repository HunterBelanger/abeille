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
#include <tallies/tallies.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <exception>
#include <sstream>
#include <string>
#include <vector>

Tallies::Tallies()
    : total_weight(0.),
      k_col_score(0.),
      k_abs_score(0.),
      k_trk_score(0.),
      leak_score(0.),
      k_tot_score(0.),
      mig_area_score(0.),
      k_col_vec(),
      k_abs_vec(),
      k_trk_vec(),
      leak_vec(),
      mig_vec(),
      k_col(1.0),
      k_col_avg(0.),
      k_col_var(0.),
      k_abs(1.),
      k_abs_avg(0.),
      k_abs_var(0.),
      k_trk(1.),
      k_trk_avg(0.),
      k_trk_var(0.),
      leak(0.0),
      leak_avg(0.0),
      leak_var(0.0),
      k_tot(1.0),
      k_tot_avg(0.),
      k_tot_var(0.),
      mig(0.),
      mig_avg(0.),
      mig_var(0.),
      collision_mesh_tallies_(),
      track_length_mesh_tallies_(),
      source_mesh_tallies_(),
      noise_source_mesh_tallies_(),
      new_itally_collision_(),
      new_itally_track_length_() {}

Tallies& Tallies::instance() {
  static Tallies tallies;

  return tallies;
}

void Tallies::allocate_batch_arrays(std::size_t nbatches) {
  k_col_vec.reserve(nbatches);
  k_abs_vec.reserve(nbatches);
  k_trk_vec.reserve(nbatches);
  leak_vec.reserve(nbatches);
  mig_vec.reserve(nbatches);
}

void Tallies::verify_track_length_tallies(bool track_length_transporter) const {
  if (track_length_transporter == false) {
    for (const auto& tally : track_length_mesh_tallies_) {
      if (tally->quantity() != MeshTally::Quantity::Flux &&
          tally->quantity() != MeshTally::Quantity::RealFlux &&
          tally->quantity() != MeshTally::Quantity::ImgFlux) {
        std::stringstream mssg;
        mssg << "Estimator for tally " << tally->name()
             << " is incompatible with the selected transport operator.";
        fatal_error(mssg.str());
      }
    }
  }
}

void Tallies::add_position_filter_set(std::size_t id,
                                      std::shared_ptr<PositionFilter> filter) {
  position_filter_set_[id] = filter;
}

void Tallies::add_cartesian_filter_set(
    std::size_t id, std::shared_ptr<CartesianFilter> filter) {
  cartesian_filter_set_[id] = filter;
}
void Tallies::add_cylinder_filter_set(
    std::size_t id, std::shared_ptr<CylinderPositionFilter> filter) {
  cylinder_position_set_[id] = filter;
}
void Tallies::add_energy_filter_set(std::size_t id,
                                    std::shared_ptr<EnergyFilter> filter) {
  energy_filter_set_[id] = filter;
}

void Tallies::add_ITally(std::shared_ptr<ITally> new_tally) {
  new_tally->set_net_weight(total_weight);

  if (new_tally->estimator() == ITally::Estimator::Collision) {
    new_itally_collision_.push_back(new_tally);
  } else if (new_tally->estimator() == ITally::Estimator::TrackLength) {
    new_itally_track_length_.push_back(new_tally);
  }
}

void Tallies::add_collision_mesh_tally(
    std::shared_ptr<CollisionMeshTally> cetally) {
  cetally->set_net_weight(total_weight);
  collision_mesh_tallies_.push_back(cetally);
}

void Tallies::add_track_length_mesh_tally(
    std::shared_ptr<TrackLengthMeshTally> tltally) {
  tltally->set_net_weight(total_weight);
  track_length_mesh_tallies_.push_back(tltally);
}

void Tallies::add_source_mesh_tally(std::shared_ptr<SourceMeshTally> stally) {
  stally->set_net_weight(total_weight);
  source_mesh_tallies_.push_back(stally);
}

void Tallies::add_noise_source_mesh_tally(
    std::shared_ptr<SourceMeshTally> stally) {
  stally->set_net_weight(total_weight);
  noise_source_mesh_tallies_.push_back(stally);
}

void Tallies::score_k_col(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  k_col_score += scr;
}

void Tallies::score_k_abs(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  k_abs_score += scr;
}

void Tallies::score_k_trk(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  k_trk_score += scr;
}

void Tallies::score_k_tot(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  k_tot_score += scr;
}

void Tallies::score_leak(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  leak_score += scr;
}

void Tallies::score_mig_area(double scr) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
  mig_area_score += scr;
}

void Tallies::clear_generation() {
  k_col_score = 0.;
  k_abs_score = 0.;
  k_trk_score = 0.;
  k_tot_score = 0.;
  leak_score = 0.;
  mig_area_score = 0.;

  for (auto& tally : collision_mesh_tallies_) tally->clear_generation();

  for (auto& tally : track_length_mesh_tallies_) tally->clear_generation();

  for (auto& tally : source_mesh_tallies_) tally->clear_generation();

  for (auto& tally : noise_source_mesh_tallies_) tally->clear_generation();

  for (auto& tally : new_itally_collision_) tally->clear_generation();

  for (auto& tally : new_itally_track_length_) tally->clear_generation();
}

void Tallies::calc_gen_values() {
  // Must first perform a reduction to get all the score
  // contributions from all processes.
  mpi::Allreduce_sum(k_col_score);
  mpi::Allreduce_sum(k_abs_score);
  mpi::Allreduce_sum(k_trk_score);
  mpi::Allreduce_sum(leak_score);
  mpi::Allreduce_sum(k_tot_score);
  mpi::Allreduce_sum(mig_area_score);

  k_col = k_col_score / total_weight;
  k_abs = k_abs_score / total_weight;
  k_trk = k_trk_score / total_weight;
  leak = leak_score / total_weight;
  k_tot = k_tot_score / total_weight;
  mig = mig_area_score / total_weight;

  k_col_vec.push_back(k_col);
  k_abs_vec.push_back(k_abs);
  k_trk_vec.push_back(k_trk);
  leak_vec.push_back(leak);
  mig_vec.push_back(mig);
}

void Tallies::record_generation(double multiplier) {
  gen++;

  update_avg_and_var(k_col, k_col_avg, k_col_var);
  update_avg_and_var(k_abs, k_abs_avg, k_abs_var);
  update_avg_and_var(k_trk, k_trk_avg, k_trk_var);
  update_avg_and_var(leak, leak_avg, leak_var);
  update_avg_and_var(k_tot, k_tot_avg, k_tot_var);
  update_avg_and_var(mig, mig_avg, mig_var);

  for (auto& tally : collision_mesh_tallies_)
    tally->record_generation(multiplier);

  for (auto& tally : track_length_mesh_tallies_)
    tally->record_generation(multiplier);

  for (auto& tally : source_mesh_tallies_) tally->record_generation(multiplier);

  for (auto& tally : noise_source_mesh_tallies_)
    tally->record_generation(multiplier);

  for (auto& new_tally : new_itally_collision_)
    new_tally->record_generation(multiplier);

  for (auto& new_tally : new_itally_track_length_)
    new_tally->record_generation(multiplier);
}

void Tallies::write_tallies(bool track_length_compatible) {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Make results group
  auto results = h5.createGroup("results");

  // Record the active number of generations
  results.createAttribute("generations", gen);

  // Write vectors to hdf5
  if (settings::sim_mode == settings::SimMode::KEFF) {
    results.createDataSet("kcol", k_col_vec);
    if (track_length_compatible) {
      // The TLE is no good if we aren't using surface tracking
      results.createDataSet("ktrk", k_trk_vec);
    }
    if (settings::energy_mode == settings::EnergyMode::CE) {
      // In MG mode, kcol = kabs, so don't wrie the same data twice.
      results.createDataSet("kabs", k_abs_vec);
    }

    if (gen > 0) {
      // Only write these if we started recording them
      results.createAttribute("kcol-avg", kcol_avg());
      results.createAttribute("kcol-std", kcol_err());

      if (track_length_compatible) {
        // The TLE is no good if we aren't using surface tracking
        results.createAttribute("ktrk-avg", ktrk_avg());
        results.createAttribute("ktrk-std", ktrk_err());
      }
      if (settings::energy_mode == settings::EnergyMode::CE) {
        // In MG mode, kcol = kabs, so don't wrie the same data twice.
        results.createAttribute("kabs-avg", kabs_avg());
        results.createAttribute("kabs-std", kabs_err());
      }
    }
  }
  results.createDataSet("leakage", leak_vec);
  results.createDataSet("mig-area", mig_vec);
  if (gen > 0) {
    results.createAttribute("leakage-avg", leakage_avg());
    results.createAttribute("leakage-str", leakage_err());

    results.createAttribute("mig-area-avg", mig_area_avg());
    results.createAttribute("mig-area-std", mig_area_err());
  }

  // Write all mesh tallies
  if (gen > 0) {
    for (auto& tally : collision_mesh_tallies_) tally->write_tally();

    for (auto& tally : track_length_mesh_tallies_) tally->write_tally();

    for (auto& tally : source_mesh_tallies_) tally->write_tally();

    for (auto& tally : noise_source_mesh_tallies_) tally->write_tally();

    for (auto& new_tally : new_itally_collision_) new_tally->write_tally();

    for (auto& new_tally : new_itally_track_length_) new_tally->write_tally();
  }
}

void Tallies::update_avg_and_var(double x, double& x_avg, double& x_var) {
  double dgen = static_cast<double>(gen);
  double x_avg_old = x_avg;
  double x_var_old = x_var;

  // Update average
  x_avg = x_avg_old + (x - x_avg_old) / (dgen);

  // Update variance
  if (gen > 1)
    x_var = x_var_old + ((x - x_avg_old) * (x - x_avg_old) / (dgen)) -
            ((x_var_old) / (dgen - 1.));
}

void make_tallies_filter(Tallies& tallies, const YAML::Node& node) {
  if (!node["tally-filters"]) {
    fatal_error("tally-filter are not given to evaluate the tallies.");
  }

  // position-filters
  if (node["tally-filters"]["position-filters"]) {
    const YAML::Node& position_filter_nodes =
        node["tally-filters"]["position-filters"];
    for (std::size_t i = 0; i < position_filter_nodes.size(); i++) {
      if (!position_filter_nodes[i]["id"] ||
          !position_filter_nodes[i]["id"].IsScalar()) {
        fatal_error(
            "Invalid \"id\" is given for position filter in tally-filters, the "
            "possible error entry is " +
            std::to_string(i) + " from top.");
      }
      std::size_t id = position_filter_nodes[i]["id"].as<std::size_t>();

      if (!position_filter_nodes[i]["position-filter-type"] ||
          !position_filter_nodes[i]["position-filter-type"].IsScalar()) {
        fatal_error("Invalid entry for the position-filter-type for " +
                    std::to_string(id) + " id.");
      }
      std::string position_filter_name =
          position_filter_nodes[i]["position-filter-type"].as<std::string>();

      if (position_filter_name == "cylinder-filter") {
        auto cylinder_filter =
            make_cylinder_position_filter(position_filter_nodes[i]);
        tallies.add_cylinder_filter_set(id, cylinder_filter);
        tallies.add_position_filter_set(id, cylinder_filter);
      } else {
        auto cartesian_filter = make_cartesian_filter(position_filter_nodes[i]);
        tallies.add_cartesian_filter_set(id, cartesian_filter);
        tallies.add_position_filter_set(id, cartesian_filter);
      }
    }
  }

  // energy-filters
  if (node["tally-filters"]["energy-filters"]) {
    const YAML::Node& energy_filter_node =
        node["tally-filters"]["energy-filters"];
    for (std::size_t i = 0; i < energy_filter_node.size(); i++) {
      if (!energy_filter_node[i]["id"] ||
          energy_filter_node[i]["id"].IsScalar()) {
        fatal_error(
            "Invalid \"id\" is given for energy filter in tally-filters, the "
            "possible error entry is " +
            std::to_string(i) + " from top.");
      }

      std::size_t id = energy_filter_node[i]["id"].as<std::size_t>();
      std::shared_ptr<EnergyFilter> filter = make_energy_filter(node);
      tallies.add_energy_filter_set(id, filter);
    }
  }
}

void add_mesh_tally(Tallies& tallies, const YAML::Node& node) {
  // First get type of estimator. Default is collision
  if (!node["type"]) {
    std::string estimator_str = "collision";
    if (node["estimator"]) {
      estimator_str = node["estimator"].as<std::string>();
    }

    if (estimator_str == "collision") {
      tallies.add_collision_mesh_tally(make_collision_mesh_tally(node));
    } else if (estimator_str == "track-length") {
      tallies.add_track_length_mesh_tally(make_track_length_mesh_tally(node));
    } else if (estimator_str == "source") {
      auto stally = make_source_mesh_tally(node);
      if (stally->noise_like_score()) {
        tallies.add_noise_source_mesh_tally(stally);
      } else {
        tallies.add_source_mesh_tally(stally);
      }
    } else {
      fatal_error("Unknown estimator type of \"" + estimator_str + "\".");
    }
  } else {
    make_itally(tallies, node);
  }
}
