/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#ifndef TALLIES_H
#define TALLIES_H

#include <yaml-cpp/yaml.h>

#include <materials/material_helper.hpp>
#include <simulation/collision_mesh_tally.hpp>
#include <simulation/particle.hpp>
#include <simulation/source_mesh_tally.hpp>
#include <simulation/track_length_mesh_tally.hpp>
#include <string>

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
