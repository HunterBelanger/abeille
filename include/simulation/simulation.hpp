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
#ifndef SIMULATION_H
#define SIMULATION_H

#include <geometry/geometry.hpp>
#include <simulation/entropy.hpp>
#include <simulation/source.hpp>
#include <simulation/tallies.hpp>
#include <simulation/transporter.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <sstream>

class Simulation {
 public:
  Simulation(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
             std::vector<std::shared_ptr<Source>> srcs);
  virtual ~Simulation() = default;

  virtual void initialize() = 0;
  virtual void run() = 0;

  virtual void premature_kill() = 0;

  // Method to sample sources
  std::vector<Particle> sample_sources(std::size_t N);

  // Methods to set entropies
  void set_p_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    p_pre_entropy = entrpy;
  }
  void set_n_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    n_pre_entropy = entrpy;
  }
  void set_t_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    t_pre_entropy = entrpy;
  }

  void set_p_post_entropy(std::shared_ptr<Entropy> entrpy) {
    p_post_entropy = entrpy;
  }
  void set_n_post_entropy(std::shared_ptr<Entropy> entrpy) {
    n_post_entropy = entrpy;
  }
  void set_t_post_entropy(std::shared_ptr<Entropy> entrpy) {
    t_post_entropy = entrpy;
  }

  bool signaled = false;
  bool terminate = false;

 protected:
  std::shared_ptr<Tallies> tallies;
  std::shared_ptr<Transporter> transporter;
  std::vector<std::shared_ptr<Source>> sources;

  Timer simulation_timer;

  uint64_t histories_counter = 0;
  uint64_t global_histories_counter = 0;
  uint64_t transported_histories = 0;

  // All entropy bins possible
  std::shared_ptr<Entropy> p_pre_entropy = nullptr;
  std::shared_ptr<Entropy> n_pre_entropy = nullptr;
  std::shared_ptr<Entropy> t_pre_entropy = nullptr;
  std::vector<double> p_pre_entropy_vec;
  std::vector<double> n_pre_entropy_vec;
  std::vector<double> t_pre_entropy_vec;

  std::shared_ptr<Entropy> p_post_entropy = nullptr;
  std::shared_ptr<Entropy> n_post_entropy = nullptr;
  std::shared_ptr<Entropy> t_post_entropy = nullptr;
  std::vector<double> p_post_entropy_vec;
  std::vector<double> n_post_entropy_vec;
  std::vector<double> t_post_entropy_vec;

  std::vector<double> empty_entropy_frac_vec;

  void sync_signaled();
  void sync_banks(std::vector<uint64_t>& nums,
                  std::vector<BankedParticle>& bank);

  void write_source(std::vector<Particle>& bank) const;

};  // Simulation

#endif  // MG_SIMULATION_H
