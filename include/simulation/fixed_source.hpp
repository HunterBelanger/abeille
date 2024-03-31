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
#ifndef FIXED_SOURCE_H
#define FIXED_SOURCE_H

#include <simulation/simulation.hpp>
#include <source/source.hpp>
#include <utils/timer.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class FixedSource : public Simulation {
 public:
  FixedSource(std::shared_ptr<IParticleMover> i_pm) : Simulation(i_pm) {}
  ~FixedSource() = default;

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;
  void write_output_info() const override final;

  void set_nparticles(std::size_t np) { nparticles = np; }
  void set_nbatches(std::size_t nb) { nbatches = nb; }
  void add_source(std::shared_ptr<Source> src);

 private:
  std::vector<std::shared_ptr<Source>> sources{};
  std::vector<double> batch_times{};
  Timer batch_timer;
  std::size_t batch = 0;
  std::size_t nparticles, nbatches;

  void check_time(std::size_t gen);
  bool out_of_time(std::size_t gen);
  void print_header();
  void batch_output(std::size_t gen);
  void mpi_setup();
  void mpi_advance();
};

std::shared_ptr<FixedSource> make_fixed_source(const YAML::Node& sim);

#endif
