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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <simulation/fixed_source.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/rng.hpp>
#include <utils/timer.hpp>

#include "utils/mpi.hpp"
#include "utils/settings.hpp"

void FixedSource::initialize() {}

void FixedSource::print_header() {
  std::shared_ptr<Output> out = Output::instance();

  // Change the header for the output based on wether or not the entropy is
  // being calculated
  out->write("\n");
  out->write("  Gen   leakage   average +/- err\n");
  out->write(" ------------------------------------\n");
  //           1120   1.23456   1.23456 +/- 0.00023
}

void FixedSource::mpi_setup() {
  // Calculate the base number of particles per node to run
  uint64_t base_particles_per_node =
      static_cast<uint64_t>(settings::nparticles / mpi::size);
  uint64_t remainder = static_cast<uint64_t>(settings::nparticles % mpi::size);

  // Set the base number of particles per node in the node_nparticles vector
  mpi::node_nparticles.resize(static_cast<std::size_t>(mpi::size),
                              base_particles_per_node);

  // Distribute the remainder particles amonst the nodes. There are at most
  // mpi::size-1 remainder particles, so we will distribute them untill we
  // have no more.
  for (std::size_t rank = 0; rank < remainder; rank++) {
    mpi::node_nparticles[rank]++;
  }

  // Now we need to make sure that the history_counter for each node is at
  // the right starting location.
  for (int lower_rank = 0; lower_rank < mpi::rank; lower_rank++) {
    histories_counter +=
        mpi::node_nparticles.at(static_cast<std::size_t>(lower_rank));
  }
}

void FixedSource::mpi_advance() {
  histories_counter +=
      static_cast<uint64_t>(settings::nparticles) -
      mpi::node_nparticles[static_cast<std::size_t>(mpi::rank)];
}

void FixedSource::run() {
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Running Fixed-Source Problem...\n");
  out->write(" NPARTICLES: " + std::to_string(settings::nparticles) + ", ");
  out->write(" NGENERATIONS: " + std::to_string(settings::ngenerations) + "\n");

  // Calculate the number of particles that each node should run.
  mpi_setup();

  // Make sure the settings is set to converged, as we dont need to
  // allow for source convergence in fixed source calculations.
  settings::converged = true;

  // Initialize vectors to hold particles
  std::vector<Particle> bank;

  // Start timer
  simulation_timer.reset();
  mpi::synchronize();
  simulation_timer.start();

  // Get the number of particles that this node should run
  uint64_t node_nparticles =
      mpi::node_nparticles[static_cast<std::size_t>(mpi::rank)];

  for (int g = 1; g <= settings::ngenerations; g++) {
    gen = g;

    // First, sample the sources and place into bank
    bank = this->sample_sources(node_nparticles);

    // Now transport all particles. In fixed-source mode, this should
    // return and empty vector !
    auto fission_bank = transporter->transport(bank);
    if (!fission_bank.empty()) {
      std::string mssg = "Returned bank not empty on fixed-source transport.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
    mpi::synchronize();

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

    // Advance history counters so that each node starts at the right location
    // for the next batch.
    mpi_advance();

    // Check if signal has been sent after generation keff has been
    // recorded, and cancellation has occured. Otherwize source file
    // will be larger than necessary, and wrong number of gens will be
    // in output file based on number averaged for tallies.
    sync_signaled();
    if (signaled) premature_kill();

    // Check if we have enough time to finish the simulation. If not,
    // stop now.
    check_time(g);
  }

  // Stop timer
  simulation_timer.stop();
  out->write("\n Total Simulation Time: " +
             std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  // Write the final results of all estimators
  out->write("\n");
  std::stringstream output;
  output << std::fixed << std::setprecision(6);
  output << " leakage = " << tallies->leakage_avg() << " +/- "
         << tallies->leakage_err() << "\n";
  out->write(output.str());

  // Write saved warnings
  out->write_saved_warnings();

  // Save other outputs
  out->write("\n");

  // Write flux file
  tallies->write_tallies();
}

void FixedSource::generation_output(int gen) {
  std::shared_ptr<Output> out = Output::instance();
  std::stringstream output;

  if (gen == 1) print_header();

  output << std::fixed << std::setw(5) << std::setfill(' ') << gen << "   "
         << std::setprecision(5) << tallies->leakage();
  output << "   " << tallies->leakage_avg() << " +/- " << tallies->leakage_err()
         << "\n";
  out->write(output.str());
}

void FixedSource::premature_kill() {
  // See if the user really wants to kill the program
  std::shared_ptr<Output> out = Output::instance();
  bool user_said_kill = false;

  if (mpi::rank == 0 && terminate == false) {
    std::string response;
    std::cout << "\n Do you really want to stop the simulation ? (y/N) => ";
    std::cin >> response;
    if (response == "y" || response == "Y") user_said_kill = true;
  }

  mpi::Bcast<bool>(user_said_kill, 0);

  if (user_said_kill || terminate) {
    out->write("\n Simulation has been stopped by user. Cleaning up...\n");
    // Looks like they really want to kill it.
    // Save everthing for them, and abort.

    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  } else if (mpi::rank == 0) {
    std::cout << "\n";
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
}

bool FixedSource::out_of_time(int gen) {
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

void FixedSource::check_time(int gen) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(gen);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  }
}
