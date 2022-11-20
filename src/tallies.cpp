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
#include <exception>
#include <simulation/tallies.hpp>
#include <string>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <vector>

Tallies::Tallies(double tot_wgt)
    : total_weight(tot_wgt),
      k_col_score(0.),
      k_abs_score(0.),
      k_trk_score(0.),
      leak_score(0.),
      k_tot_score(0.),
      mig_area_score(0.),
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
      noise_source_mesh_tallies_() {}

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
#ifdef _OPENMP
#pragma omp atomic
#endif
  k_col_score += scr;
}

void Tallies::score_k_abs(double scr) {
#ifdef _OPENMP
#pragma omp atomic
#endif
  k_abs_score += scr;
}

void Tallies::score_k_trk(double scr) {
#ifdef _OPENMP
#pragma omp atomic
#endif
  k_trk_score += scr;
}

void Tallies::score_k_tot(double scr) {
#ifdef _OPENMP
#pragma omp atomic
#endif
  k_tot_score += scr;
}

void Tallies::score_leak(double scr) {
#ifdef _OPENMP
#pragma omp atomic
#endif
  leak_score += scr;
}

void Tallies::score_mig_area(double scr) {
#ifdef _OPENMP
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
}

void Tallies::write_tallies() {
  for (auto& tally : collision_mesh_tallies_) tally->write_tally();

  for (auto& tally : track_length_mesh_tallies_) tally->write_tally();

  for (auto& tally : source_mesh_tallies_) tally->write_tally();

  for (auto& tally : noise_source_mesh_tallies_) tally->write_tally();
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

void add_mesh_tally(Tallies& tallies, const YAML::Node& node) {
  // First get type of estimator. Default is collision
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
    std::string mssg = "Unknown estimator type of \"" + estimator_str + "\".";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}
