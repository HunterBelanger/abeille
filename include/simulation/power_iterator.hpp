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
#ifndef POWER_ITERATOR_H
#define POWER_ITERATOR_H

#include <cancelator/cancelator.hpp>
#include <simulation/entropy.hpp>
#include <simulation/simulation.hpp>
#include <source/source.hpp>
#include <utils/timer.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

class PowerIterator : public Simulation {
 public:
  PowerIterator(std::shared_ptr<IParticleMover> i_pm) : Simulation(i_pm) {};
  ~PowerIterator() = default;

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;
  void write_output_info() const override final;

  void set_cancelator(std::shared_ptr<Cancelator> cncl);
  void set_entropy(const Entropy& entropy);
  void set_combing(bool comb) { combing = comb; }
  void set_families(bool fams) { calc_families = fams; }
  void set_pair_distance(bool prdist) { pair_distance_sqrd = prdist; }
  void set_empty_entropy_bins(bool eeb) { empty_entropy_bins = eeb; }
  void set_save_source(bool save_src) { save_source = save_src; }
  void set_in_source_file(const std::string& sf) { in_source_file_name = sf; }
  void add_source(std::shared_ptr<Source> src);
  void set_nparticles(std::size_t np) { nparticles = np; }
  void set_ngenerations(std::size_t ng) { ngenerations = ng; }
  void set_nignored(std::size_t ni) { nignored = ni; }

 private:
  std::vector<std::shared_ptr<Source>> sources{};
  std::vector<Particle> bank{};
  std::vector<std::size_t> families_vec{};
  std::vector<int> Nnet_vec{}, Npos_vec{}, Nneg_vec{}, Ntot_vec{};
  std::vector<int> Wnet_vec{}, Wpos_vec{}, Wneg_vec{}, Wtot_vec{};
  std::vector<double> r_sqrd_vec{};
  std::vector<double> t_pre_entropy_vec{}, p_pre_entropy_vec{},
      n_pre_entropy_vec{};
  std::vector<double> t_post_entropy_vec{}, p_post_entropy_vec{},
      n_post_entropy_vec{};
  std::vector<double> empty_entropy_frac_vec{};
  std::string in_source_file_name{};
  std::unordered_set<uint64_t> families_set{};
  std::vector<double> generation_times{};
  Timer generation_timer;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  std::shared_ptr<Entropy> t_pre_entropy = nullptr;
  std::shared_ptr<Entropy> p_pre_entropy = nullptr;
  std::shared_ptr<Entropy> n_pre_entropy = nullptr;
  std::shared_ptr<Entropy> t_post_entropy = nullptr;
  std::shared_ptr<Entropy> p_post_entropy = nullptr;
  std::shared_ptr<Entropy> n_post_entropy = nullptr;
  std::size_t ngenerations, nignored, nparticles;
  std::size_t gen = 0;
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;
  double r_sqrd = 0.;
  bool converged = false;
  bool combing = false;
  bool calc_families = false;
  bool pair_distance_sqrd = false;
  bool empty_entropy_bins = false;
  bool save_source = false;

  void check_time(std::size_t gen);

  bool out_of_time(std::size_t gen);

  void generation_output();

  void write_entropy_families_etc_to_results() const;

  void save_weights();

  // Entropy calculation
  void compute_pre_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void compute_post_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void zero_entropy();

  // Pair-distance sqrd calculation
  void compute_pair_dist_sqrd(const std::vector<BankedParticle>& next_gen);

  void print_header();

  void sample_source_from_sources();
  void load_source_from_file();

};  // PowerIterator

std::shared_ptr<PowerIterator> make_power_iterator(const YAML::Node& sim);

#endif  // MG_POWER_ITERATOR_H
