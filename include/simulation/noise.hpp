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
#ifndef NOISE_H
#define NOISE_H

#include <cancelator/cancelator.hpp>
#include <noise_source/noise_maker.hpp>
#include <simulation/simulation.hpp>

#include <memory>

class Noise : public Simulation {
 public:
  Noise(std::shared_ptr<IParticleMover> i_pm,
        std::vector<std::shared_ptr<Source>> srcs, NoiseMaker noise_mkr)
      : Simulation(i_pm, srcs),
        noise_maker(noise_mkr),
        bank(),
        noise_bank(),
        noise_timer(),
        cancellation_timer(),
        power_iteration_timer(),
        convergence_timer(),
        noise_batch_timer() {}

  Noise(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
        std::vector<std::shared_ptr<Source>> srcs,
        std::shared_ptr<Cancelator> cancel, NoiseMaker noise_mkr)
      : Simulation(i_t, i_tr, srcs),
        noise_maker(noise_mkr),
        bank(),
        noise_bank(),
        noise_timer(),
        cancellation_timer(),
        power_iteration_timer(),
        convergence_timer(),
        noise_batch_timer(),
        cancelator(cancel) {}

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;

 private:
  NoiseMaker noise_maker;
  std::vector<Particle> bank;
  std::vector<BankedParticle> noise_bank;
  Timer noise_timer;
  Timer cancellation_timer;
  Timer power_iteration_timer;
  Timer convergence_timer;
  Timer noise_batch_timer;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  int noise_batch = 0;  // Counter for number of noise batches
  int pi_gen = 0;       // Counter for number of power iteration generations
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;

  // Private helper methods
  void compute_pre_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void compute_post_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void zero_entropy();

  bool out_of_time(int gen);
  void check_time(int gen);

  void normalize_weights(std::vector<BankedParticle>& next_gen);
  void perform_regional_cancellation(std::vector<BankedParticle>& bank);

  void power_iteration(bool sample_noise);
  void pi_generation_output();

  void noise_simulation();
  void noise_output();

  void print_header() const;

  void sample_source_from_sources();
  void load_source_from_file();
};

#endif
