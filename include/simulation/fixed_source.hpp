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
#include <utils/timer.hpp>

#include <vector>

class FixedSource : public Simulation {
 public:
  FixedSource(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
              std::vector<std::shared_ptr<Source>> src)
      : Simulation(i_t, i_tr, src) {}
  ~FixedSource() = default;

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;

 private:
  std::vector<double> batch_times{};
  Timer batch_timer;
  int gen = 0;
  void check_time(int gen);
  bool out_of_time(int gen);
  void print_header();
  void generation_output(int gen);
  void mpi_setup();
  void mpi_advance();
};

#endif
