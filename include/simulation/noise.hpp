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
#include <simulation/entropy.hpp>
#include <simulation/simulation.hpp>
#include <source/source.hpp>
#include <utils/noise_parameters.hpp>

#include <memory>
#include <numeric>

class Noise : public Simulation {
 public:
  Noise(std::shared_ptr<INoiseParticleMover> i_pm,
        std::shared_ptr<IParticleMover> i_pipm, NoiseParameters noise_parms,
        NoiseMaker noise_mkr)
      : Simulation(i_pipm),
        noise_params(noise_parms),
        noise_maker(noise_mkr),
        noise_particle_mover(i_pm) {}

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;
  void write_output_info() const override final;

  void set_in_source_file(const std::string& sf) { in_source_file_name = sf; }
  void add_source(std::shared_ptr<Source> src);
  void set_nparticles(std::size_t np) { nparticles = np; }
  void set_nbatches(std::size_t nb) { nbatches = nb; }
  void set_nignored(std::size_t ni) { nignored = ni; }
  void set_nskip(std::size_t ns) { nskip = ns; }
  void set_ncancel_noise_gens(std::size_t ncng) { ncancel_noise_gens = ncng; }
  void set_normalize_noise_source(bool nns) { normalize_noise_source_ = nns; }
  void set_entropy(const Entropy& entropy);
  void set_cancelator(std::shared_ptr<Cancelator> cncl);
  void set_combing(bool comb) { combing = comb; }
  void set_regional_cancellation(bool rc);
  void set_regional_cancellation_noise(bool rcn);

 private:
  NoiseParameters noise_params;
  NoiseMaker noise_maker;
  std::vector<Particle> bank{};
  std::vector<BankedParticle> noise_bank{};
  std::vector<std::shared_ptr<Source>> sources{};
  std::string in_source_file_name{};
  Timer noise_timer{};
  Timer cancellation_timer{};
  Timer power_iteration_timer{};
  Timer convergence_timer{};
  Timer noise_batch_timer{};
  std::shared_ptr<INoiseParticleMover> noise_particle_mover = nullptr;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  std::shared_ptr<Entropy> t_pre_entropy = nullptr;
  std::shared_ptr<Entropy> p_pre_entropy = nullptr;
  std::shared_ptr<Entropy> n_pre_entropy = nullptr;
  std::shared_ptr<Entropy> t_post_entropy = nullptr;
  std::shared_ptr<Entropy> p_post_entropy = nullptr;
  std::shared_ptr<Entropy> n_post_entropy = nullptr;
  std::size_t noise_batch = 0;  // Counter for number of noise batches
  std::size_t pi_gen = 0;  // Counter for number of power iteration generations
  std::size_t nbatches = 0;
  std::size_t nparticles = 0;
  std::size_t nignored = 0;
  std::size_t nskip = 3;
  std::size_t ncancel_noise_gens = std::numeric_limits<std::size_t>::max();
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;
  bool converged = false;
  bool combing = false;
  bool regional_cancellation_ = false;
  bool regional_cancellation_noise_ = false;
  bool normalize_noise_source_ = true;

  // Private helper methods
  void compute_pre_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void compute_post_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void zero_entropy();

  bool out_of_time(std::size_t batch);
  void check_time(std::size_t batch);

  void power_iteration(bool sample_noise);
  void pi_generation_output();

  void noise_simulation();
  void noise_output();

  void print_header() const;

  void sample_source_from_sources();
  void load_source_from_file();
};

std::shared_ptr<Noise> make_noise_simulator(const YAML::Node& sim);

#endif
