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

#include <simulation/cancelator.hpp>
#include <simulation/simulation.hpp>

#include <unordered_set>
#include <vector>

class PowerIterator : public Simulation {
 public:
  PowerIterator(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
                std::vector<std::shared_ptr<Source>> src)
      : Simulation(i_t, i_tr, src), bank(), r_sqrd_vec(){};
  PowerIterator(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
                std::vector<std::shared_ptr<Source>> src,
                std::shared_ptr<Cancelator> cncl)
      : Simulation(i_t, i_tr, src), bank(), cancelator(cncl), r_sqrd_vec(){};
  ~PowerIterator() = default;

  void initialize() override final;

  void run() override final;

  void premature_kill() override final;

 private:
  std::vector<Particle> bank;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  std::unordered_set<uint64_t> families = {};
  std::vector<std::size_t> families_vec = {};
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;
  std::vector<int> Nnet_vec{}, Npos_vec{}, Nneg_vec{}, Ntot_vec{};
  std::vector<int> Wnet_vec{}, Wpos_vec{}, Wneg_vec{}, Wtot_vec{};
  int gen = 0;
  double r_sqrd = 0.;
  std::vector<double> r_sqrd_vec;

  void check_time(int gen);

  bool out_of_time(int gen);

  void generation_output();

  void write_entropy_families_etc_to_results() const;

  void normalize_weights(std::vector<BankedParticle>& next_gen);

  void perform_regional_cancellation(std::vector<BankedParticle>& next_gen);

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

#endif  // MG_POWER_ITERATOR_H
