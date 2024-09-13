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

#include <cancelator/cancelator.hpp>
#include <simulation/particle.hpp>
#include <simulation/particle_mover.hpp>
#include <source/source.hpp>
#include <utils/timer.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <sstream>
#include <vector>

class Simulation {
 public:
  Simulation(std::shared_ptr<IParticleMover> i_pm);
  virtual ~Simulation() = default;

  virtual void initialize() = 0;
  virtual void run() = 0;
  virtual void premature_kill() = 0;
  virtual void write_output_info() const = 0;

  std::vector<Particle> sample_sources(
      const std::vector<std::shared_ptr<Source>>& sources, std::size_t N);

  bool signaled = false;
  bool terminate = false;

 protected:
  std::shared_ptr<IParticleMover> particle_mover;

  Timer simulation_timer;

  uint64_t histories_counter = 0;
  uint64_t global_histories_counter = 0;
  uint64_t transported_histories = 0;

  void sync_signaled();
  void sync_banks(std::vector<uint64_t>& nums,
                  std::vector<BankedParticle>& bank);

  void normalize_weights(std::size_t nparticles,
                         std::vector<BankedParticle>& next_gen, int& Npos,
                         int& Nneg, int& Nnet, int& Ntot, int& Wpos, int& Wneg,
                         int& Wnet, int& Wtot);

  void comb_particles(std::vector<BankedParticle>& next_gen, int& Npos,
                      int& Nneg, int& Nnet, int& Ntot);

  void perform_regional_cancellation(std::shared_ptr<Cancelator>& cancelator,
                                     std::vector<BankedParticle>& next_gen);

  void write_source(std::vector<Particle>& bank) const;

};  // Simulation

std::shared_ptr<Simulation> make_simulation(const YAML::Node& input);

#endif  // MG_SIMULATION_H
