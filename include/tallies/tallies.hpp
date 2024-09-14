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
#ifndef TALLIES_H
#define TALLIES_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <simulation/tracker.hpp>
#include <tallies/cartesian_filter.hpp>
#include <tallies/cylinder_filter.hpp>
#include <tallies/energy_filter.hpp>
#include <tallies/itally.hpp>
#include <tallies/position_filter.hpp>

#include <yaml-cpp/yaml.h>

#include <map>
#include <set>
#include <string>

struct ThreadLocalScores {
  double k_col_score = 0.;
  double k_abs_score = 0.;
  double k_trk_score = 0.;
  double k_tot_score = 0.;
  double leakage_score = 0.;
  double mig_score = 0.;

  void set_zero() {
    k_col_score = 0.;
    k_abs_score = 0.;
    k_trk_score = 0.;
    k_tot_score = 0.;
    leakage_score = 0.;
    mig_score = 0.;
  }
};

// This is a singleton which holds all the tallies for the simulation
class Tallies {
 public:
  Tallies(const Tallies&) = delete;
  Tallies(Tallies&&) = delete;
  Tallies& operator=(const Tallies&) = delete;
  Tallies& operator=(Tallies&&) = delete;
  ~Tallies() = default;

  static Tallies& instance();

  //===============================================
  // addding filters instances into the Map
  void add_position_filter(std::size_t id,
                           std::shared_ptr<PositionFilter> filter);
  void add_cartesian_filter(std::size_t id,
                            std::shared_ptr<CartesianFilter> filter);
  void add_cylinder_filter(std::size_t id,
                           std::shared_ptr<CylinderFilter> filter);
  void add_energy_filter(std::size_t id, std::shared_ptr<EnergyFilter> filter);

  //===============================================
  // Method to get the filters from Maps
  std::shared_ptr<PositionFilter> get_position_filter(std::size_t id) {
    if (position_filters_.find(id) == position_filters_.end()) {
      return nullptr;
    }
    return position_filters_[id];
  }
  std::shared_ptr<CartesianFilter> get_cartesian_filter(std::size_t id) {
    if (cartesian_filters_.find(id) == cartesian_filters_.end()) {
      return nullptr;
    }
    return cartesian_filters_[id];
  }
  std::shared_ptr<CylinderFilter> get_cylinder_filter(std::size_t id) {
    if (cylinder_filters_.find(id) == cylinder_filters_.end()) {
      return nullptr;
    }
    return cylinder_filters_[id];
  }
  std::shared_ptr<EnergyFilter> get_energy_filter(std::size_t id) {
    if (energy_filters_.find(id) == energy_filters_.end()) {
      return nullptr;
    }
    return energy_filters_[id];
  }
  void add_tally(std::shared_ptr<ITally> tally);

  void allocate_batch_arrays(std::size_t nbatches);

  void verify_track_length_tallies(bool track_length_transporter) const;

  void score_collision(const Particle& p, const Tracker& tktr,
                       MaterialHelper& mat) {
    if (scoring_ && !new_itally_collision_.empty()) {
      for (auto& tally : new_itally_collision_) {
        tally->score_collision(p, tktr, mat);
      }
    }
  }

  void score_flight(const Particle& p, const Tracker& trkr, double d_flight,
                    MaterialHelper& mat) {
    if (scoring_ && !new_itally_track_length_.empty()) {
      for (auto& tally : new_itally_track_length_) {
        tally->score_flight(p, trkr, d_flight, mat);
      }
    }
  }

  void score_source(const BankedParticle& p) {
    if (scoring_ && !new_itally_source_.empty()) {
      for (auto& tally : new_itally_source_) tally->score_source(p);
    }
  }

  void score_source(const std::vector<BankedParticle>& vp) {
    if (scoring_ && !new_itally_source_.empty()) {
      for (const auto& p : vp) {
        for (auto& tally : new_itally_source_) tally->score_source(p);
      }
    }
  }

  void score_noise_source(const BankedParticle& p) {
    if (scoring_) {
      for (auto& tally : new_itally_source_) {
        auto typ = tally->quantity().type;
        if (typ == Quantity::Type::RealSource ||
            typ == Quantity::Type::ImagSource)
          tally->score_source(p);
      }
    }
  }

  void score_noise_source(const std::vector<BankedParticle>& vp) {
    if (scoring_ && !new_itally_source_.empty()) {
      for (const auto& p : vp) {
        for (auto& tally : new_itally_source_) {
          auto typ = tally->quantity().type;
          if (typ == Quantity::Type::RealSource ||
              typ == Quantity::Type::ImagSource)
            tally->score_source(p);
        }
      }
    }
  }

  void score_k_col(double scr);
  void score_k_abs(double scr);
  void score_k_trk(double scr);
  void score_leak(double scr);
  void score_k_tot(double scr);
  void score_mig_area(double scr);

  void clear_generation();

  void calc_gen_values();

  void record_generation(double multiplier = 1.);

  double kcol() const { return k_col; }
  void set_kcol(double k) { k_col = k; }  // Only used for noise !
  double kcol_avg() const { return k_col_avg; }
  double kcol_err() const {
    return std::sqrt(k_col_var / static_cast<double>(gen));
  }

  double kabs() const { return k_abs; }
  double kabs_avg() const { return k_abs_avg; }
  double kabs_err() const {
    return std::sqrt(k_abs_var / static_cast<double>(gen));
  }

  double ktrk() const { return k_trk; }
  double ktrk_avg() const { return k_trk_avg; }
  double ktrk_err() const {
    return std::sqrt(k_trk_var / static_cast<double>(gen));
  }

  double leakage() const { return leak; }
  double leakage_avg() const { return leak_avg; }
  double leakage_err() const {
    return std::sqrt(leak_var / static_cast<double>(gen));
  }

  double ktot() const { return k_tot; }
  double ktot_avg() const { return k_tot_avg; }
  double ktot_err() const {
    return std::sqrt(k_tot_var / static_cast<double>(gen));
  }

  double alpha() const { return alpha_; }
  void set_alpha(double a) { alpha_ = a; }
  double alpha_avg() const { return alpha_avg_; }
  double alpha_err() const {
    return std::sqrt(alpha_var_ / static_cast<double>(gen));
  }

  double k_alpha() const { return k_alpha_; }
  void set_k_alpha(double a) { k_alpha_ = a; }
  double k_alpha_avg() const { return k_alpha_avg_; }
  double k_alpha_err() const {
    return std::sqrt(k_alpha_var_ / static_cast<double>(gen));
  }

  double mig_area() const { return mig; }
  double mig_area_avg() const { return mig_avg; }
  double mig_area_err() const {
    return std::sqrt(mig_var / static_cast<double>(gen));
  }

  void write_tallies(bool track_length_compatible);

  void set_total_weight(double tot_wgt) {
    total_weight = tot_wgt;

    for (auto& t : new_itally_collision_) t->set_net_weight(total_weight);
    for (auto& t : new_itally_track_length_) t->set_net_weight(total_weight);
    for (auto& t : new_itally_source_) t->set_net_weight(total_weight);
  }

  int generations() const { return gen; }

  bool scoring() const { return scoring_; }

  void set_scoring(bool scr) { scoring_ = scr; }

 private:
  Tallies();
  bool scoring_ = true;
  int gen = 0;

  double total_weight;

  double k_col_score;
  double k_abs_score;
  double k_trk_score;
  double leak_score;
  double k_tot_score;  // Weird keff for NWDT
  double mig_area_score;

  std::vector<double> k_col_vec;
  std::vector<double> k_abs_vec;
  std::vector<double> k_trk_vec;
  std::vector<double> alpha_vec;
  std::vector<double> k_alpha_vec;
  std::vector<double> leak_vec;
  std::vector<double> mig_vec;

  double k_col, k_col_avg, k_col_var;
  double k_abs, k_abs_avg, k_abs_var;
  double k_trk, k_trk_avg, k_trk_var;
  double alpha_, alpha_avg_, alpha_var_;
  double k_alpha_, k_alpha_avg_, k_alpha_var_;
  double leak, leak_avg, leak_var;
  double k_tot, k_tot_avg, k_tot_var;  // Weird keff for NWDT
  double mig, mig_avg, mig_var;

  std::vector<std::shared_ptr<ITally>>
      new_itally_collision_;  // New "Itally" vector for collision
  std::vector<std::shared_ptr<ITally>>
      new_itally_track_length_;  // New "Itally" vector for track-length
  std::vector<std::shared_ptr<ITally>>
      new_itally_source_;  // New "Itally" vector for source

  // mapes for the position, cartesian, cylinder-position, and energy-filter
  std::map<std::size_t, std::shared_ptr<PositionFilter>> position_filters_;
  std::map<std::size_t, std::shared_ptr<CartesianFilter>> cartesian_filters_;
  std::map<std::size_t, std::shared_ptr<CylinderFilter>> cylinder_filters_;
  std::map<std::size_t, std::shared_ptr<EnergyFilter>> energy_filters_;

  std::set<std::string> taken_tally_names_;

  void update_avg_and_var(double x, double& x_avg, double& x_var);

  void write_filters();

};  // Tallies

void add_tally(Tallies& tallies, const YAML::Node& node);

// for making the tallies filter
void make_tally_filters(Tallies& tallies, const YAML::Node& node);

#endif  // MG_TALLIES_H
