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
#include <simulation/modified_fixed_source.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <iomanip>

void ModifiedFixedSource::initialize() {}

void ModifiedFixedSource::print_header() {
  Output& out = Output::instance();

  // Change the header for the output based on wether or not the entropy is
  // being calculated
  out.write("\n");
  out.write("  Gen   leakage   average +/- err\n");
  out.write(" ------------------------------------\n");
  //           1120   1.23456   1.23456 +/- 0.00023
}

void ModifiedFixedSource::generation_output(int gen) {
  std::stringstream output;

  if (gen == 0) print_header();

  output << std::fixed << std::setw(5) << std::setfill(' ') << gen << "   "
         << std::setprecision(5) << tallies->leakage();
  output << "   " << tallies->leakage_avg() << " +/- " << tallies->leakage_err()
         << "\n";
  Output::instance().write(output.str());
}

void ModifiedFixedSource::run() {
  Output& out = Output::instance();
  out.write(" Running Modified-Fixed-Source Problem...\n");
  out.write(" NPARTICLES: " + std::to_string(settings::nparticles) + ", ");
  out.write(" NGENERATIONS: " + std::to_string(settings::ngenerations) + "\n");

  // Make sure the settings is set to converged, as we dont need to
  // allow for source convergence in fixed source calculations.
  settings::converged = true;

  // Initialize vectors to hold particles
  std::vector<Particle> bank;

  // Start timer
  simulation_timer.reset();
  simulation_timer.start();

  for (int g = 1; g <= settings::ngenerations; g++) {
    gen = g;

    // First, sample the sources and place into bank
    bank = this->sample_sources(static_cast<std::size_t>(settings::nparticles));

    while (!bank.empty()) {
      auto fission_bank = particle_mover->transport(bank);
      transported_histories += bank.size();

      // Cancellation may be performed on fission_bank here

      bank.clear();
      bank.reserve(fission_bank.size());
      for (auto& p : fission_bank) {
        bank.emplace_back(p.r, p.u, p.E, p.wgt, p.wgt2, histories_counter++);
        bank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
      }
    }

    // Get new values
    tallies->calc_gen_values();

    // Keep values
    tallies->record_generation();

    // Zero tallies for next generation
    tallies->clear_generation();

    // Clear particle banks
    bank.clear();

    // Output
    generation_output(g);

    // Check if signal has been sent after generation keff has been
    // recorded, and cancellation has occured. Otherwize source file
    // will be larger than necessary, and wrong number of gens will be
    // in output file based on number averaged for tallies.
    if (signaled) premature_kill();

    // Check if we have enough time to finish the simulation. If not,
    // stop now.
    check_time(g);
  }

  // Stop timer
  simulation_timer.stop();
  out.write("\n Total Simulation Time: " +
            std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  // Write the final results of all estimators
  out.write("\n");
  std::stringstream output;
  output << std::fixed << std::setprecision(6);
  output << " leakage = " << tallies->leakage_avg() << " +/- "
         << tallies->leakage_err() << "\n";
  out.write(output.str());

  // Write saved warnings
  out.write_saved_warnings();

  // Save other outputs
  out.write("\n");

  // Write flux file
  tallies->write_tallies();
}

void ModifiedFixedSource::premature_kill() {
  // See if the user really wants to kill the program
  std::string response;
  std::cout << "\n Do you really want to stop the simulation ? (y/N) => ";
  std::cin >> response;

  if (response == "y" || response == "Y") {
    Output::instance().write(
        "\n Simulation has been stopped by user. Cleaning up...\n");
    // Looks like they really want to kill it.
    // Save everthing for them, and abort.

    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

bool ModifiedFixedSource::out_of_time(int gen) {
  // Get the average time per batch
  double T_avg = simulation_timer.elapsed_time() / static_cast<double>(gen);

  // See how much time we have used so far.
  double T_used = settings::alpha_omega_timer.elapsed_time();

  // See how much time we have left
  double T_remaining = settings::max_time - T_used;

  // If the remaining time is less than 2*T_avg, than we just kill it
  // it here, so that we are sure we don't run over.
  if (T_remaining < 2. * T_avg) {
    return true;
  }

  return false;
}

void ModifiedFixedSource::check_time(int gen) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(gen);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  }
}
