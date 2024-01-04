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
#include <tallies/collision_mesh_tally.hpp>
#include <tallies/source_mesh_tally.hpp>
#include <tallies/track_length_mesh_tally.hpp>

#include <yaml-cpp/yaml.h>

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

class Tallies {
 public:
  Tallies(double tot_wgt);
  ~Tallies() = default;

  void add_collision_mesh_tally(std::shared_ptr<CollisionMeshTally> cetally);
  void add_track_length_mesh_tally(
      std::shared_ptr<TrackLengthMeshTally> tltally);
  void add_source_mesh_tally(std::shared_ptr<SourceMeshTally> stally);
  void add_noise_source_mesh_tally(std::shared_ptr<SourceMeshTally> stally);

  void score_collision(const Particle& p, MaterialHelper& mat, bool converged) {
    // Only do spacial tallies if converged
    if (converged && !collision_mesh_tallies_.empty()) {
      for (auto& tally : collision_mesh_tallies_)
        tally->score_collision(p, mat);
    }
  }

  void score_flight(const Particle& p, double d, MaterialHelper& mat,
                    bool converged) {
    if (converged && !track_length_mesh_tallies_.empty()) {
      for (auto& tally : track_length_mesh_tallies_)
        tally->score_flight(p, d, mat);
    }
  }

  void score_source(const BankedParticle& p, bool converged) {
    if (converged && !source_mesh_tallies_.empty()) {
      for (auto& tally : source_mesh_tallies_) tally->score_source(p);
    }
  }

  void score_source(const std::vector<BankedParticle>& vp, bool converged) {
    if (converged && !source_mesh_tallies_.empty()) {
      for (const auto& p : vp) {
        for (auto& tally : source_mesh_tallies_) tally->score_source(p);
      }
    }
  }

  void score_noise_source(const BankedParticle& p, bool converged) {
    if (converged && !noise_source_mesh_tallies_.empty()) {
      for (auto& tally : noise_source_mesh_tallies_) tally->score_source(p);
    }
  }

  void score_noise_source(const std::vector<BankedParticle>& vp,
                          bool converged) {
    if (converged && !noise_source_mesh_tallies_.empty()) {
      for (const auto& p : vp) {
        for (auto& tally : noise_source_mesh_tallies_) tally->score_source(p);
      }
    }
  }

  void set_keff(double k) { keff_ = k; }
  double keff() const { return keff_; }

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

  double mig_area() const { return mig; }
  double mig_area_avg() const { return mig_avg; }
  double mig_area_err() const {
    return std::sqrt(mig_var / static_cast<double>(gen));
  }

  void write_tallies();

  void set_total_weight(double tot_wgt) { total_weight = tot_wgt; }

  int generations() const { return gen; }

 private:
  int gen = 0;
  double keff_ = 1.;

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
  std::vector<double> leak_vec;
  std::vector<double> mig_vec;

  double k_col, k_col_avg, k_col_var;
  double k_abs, k_abs_avg, k_abs_var;
  double k_trk, k_trk_avg, k_trk_var;
  double leak, leak_avg, leak_var;
  double k_tot, k_tot_avg, k_tot_var;  // Weird keff for NWDT
  double mig, mig_avg, mig_var;

  std::vector<std::shared_ptr<CollisionMeshTally>> collision_mesh_tallies_;
  std::vector<std::shared_ptr<TrackLengthMeshTally>> track_length_mesh_tallies_;
  std::vector<std::shared_ptr<SourceMeshTally>> source_mesh_tallies_;
  std::vector<std::shared_ptr<SourceMeshTally>> noise_source_mesh_tallies_;

  void update_avg_and_var(double x, double& x_avg, double& x_var);

};  // Tallies

void add_mesh_tally(Tallies& tallies, const YAML::Node& node);

#endif  // MG_TALLIES_H
