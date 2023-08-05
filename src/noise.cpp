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
#include <simulation/noise.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <iomanip>
#include <numeric>
#include <sstream>

// Here, we sample a particle source exactly as we do for power iteration.
void Noise::initialize() {
  // If not using source file
  if (settings::load_source_file == false) {
    sample_source_from_sources();
  } else {
    // Load source file
    load_source_from_file();
  }

  mpi::node_nparticles_noise.resize(static_cast<std::size_t>(mpi::size), 0);
}

void Noise::load_source_from_file() {
  // Load source file
  auto h5s = H5::File(settings::in_source_file_name, H5::File::ReadOnly);

  // Get the dataset and dimensions
  auto source_ds = h5s.getDataSet("source");
  std::vector<std::size_t> dimensions = source_ds.getSpace().getDimensions();

  // Check dimensions
  if (dimensions.size() != 2 || dimensions[1] != 9) {
    fatal_error("Invalid source from file dimensions.");
  }

  // Read in array
  NDArray<double> source = NDArray<double>(dimensions);
  source_ds.read<double>(&source[0]);

  // Get number of particles
  std::size_t Nprt = source.shape()[0];

  // Calculate the base number of particles per node to run
  uint64_t base_particles_per_node =
      static_cast<uint64_t>(Nprt / static_cast<std::size_t>(mpi::size));
  uint64_t remainder =
      static_cast<uint64_t>(Nprt % static_cast<std::size_t>(mpi::size));

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

  // Each node starts reading the input source data at their histories_counter
  // index. It then reads its assigned number of particles.
  uint64_t file_start_loc = histories_counter;
  double tot_wgt = 0.;
  for (std::size_t i = 0;
       i < mpi::node_nparticles[static_cast<std::size_t>(mpi::rank)]; i++) {
    double x = source[(file_start_loc + i) * 9 + 0];
    double y = source[(file_start_loc + i) * 9 + 1];
    double z = source[(file_start_loc + i) * 9 + 2];
    double ux = source[(file_start_loc + i) * 9 + 3];
    double uy = source[(file_start_loc + i) * 9 + 4];
    double uz = source[(file_start_loc + i) * 9 + 5];
    double E = source[(file_start_loc + i) * 9 + 6];
    double w = source[(file_start_loc + i) * 9 + 7];
    double w2 = source[(file_start_loc + i) * 9 + 8];

    bank.push_back({{x, y, z}, {ux, uy, uz}, E, w, histories_counter++});
    bank.back().set_weight2(w2);
    bank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
    tot_wgt += w;
  }
  global_histories_counter = Nprt;

  mpi::Allreduce_sum(tot_wgt);

  Output::instance().write(
      " Total Weight of System: " + std::to_string(std::round(tot_wgt)) + "\n");
  settings::nparticles = static_cast<int>(std::round(tot_wgt));
  tallies->set_total_weight(std::round(tot_wgt));
}

void Noise::sample_source_from_sources() {
  Output::instance().write(" Generating source particles...\n");
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

  // Go sample the particles for this node
  bank =
      sample_sources(mpi::node_nparticles[static_cast<std::size_t>(mpi::rank)]);

  global_histories_counter += static_cast<uint64_t>(std::accumulate(
      mpi::node_nparticles.begin(), mpi::node_nparticles.end(), 0));
}

void Noise::print_header() const {
  std::stringstream output;
  output << "\n";

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(settings::ngenerations).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << "Gen"
         << "      k     ";

  if (!settings::regional_cancellation && t_pre_entropy) {
    output << "    Entropy  ";
  }

  // Add line underneath
  output << "\n " << std::string(output.str().size() - 3, '-') << "\n";

  Output::instance().write(output.str());
}

bool Noise::out_of_time(int gen) {
  // Get the average time per batch
  double T_avg = noise_batch_timer.elapsed_time() / static_cast<double>(gen);

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

void Noise::check_time(int gen) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(gen);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  }
}

void Noise::run() {
  Output& out = Output::instance();
  out.write(" Running Noise Problem...\n");
  out.write(" NPARTICLES: " + std::to_string(settings::nparticles) + ", ");
  out.write(" NBATCHES: " + std::to_string(settings::ngenerations) + ",");
  out.write(" NIGNORED: " + std::to_string(settings::nignored) + ",");
  out.write(" NSKIP: " + std::to_string(settings::nskip) + "\n\n");
  print_header();

  // Zero all entropy bins
  zero_entropy();

  // Start timer
  simulation_timer.reset();
  simulation_timer.start();

  //============================================================================
  // Run standard power iteration untill the source has converged.
  // Make sure converged is false so we don't collect statistics.
  convergence_timer.reset();
  convergence_timer.start();
  settings::converged = false;
  for (int g = 1; g <= settings::nignored; g++) {
    pi_gen = g;
    power_iteration(false);
  }

  // We now consider the source to be converged.
  settings::converged = true;

  // Make sure timers is cleared
  convergence_timer.stop();
  noise_timer.reset();
  cancellation_timer.reset();
  power_iteration_timer.reset();

  //============================================================================
  // Now we run the desired number of noise batches, which is just ngenerations
  noise_batch_timer.reset();
  noise_batch_timer.start();
  for (noise_batch = 1; noise_batch <= settings::ngenerations; noise_batch++) {
    // We first run the decorrelation PI batches, where no source is sampled.
    // This should be nskip - 1.
    // Make sure converged is false so we don't collect statistics on the
    // power iteration parts.
    power_iteration_timer.start();
    for (int skip = 0; skip < settings::nskip - 1; skip++) {
      pi_gen++;
      power_iteration(false);
    }

    // We now run the last remaining power iteration, but here we sample the
    // noise source
    pi_gen++;
    power_iteration(true);
    power_iteration_timer.stop();

    // The noise_bank is now populated with noise particles. We may now simulate
    // these noise particles.
    noise_simulation();

    if (signaled) premature_kill();

    // Check if we have enough time to finish the simulation. If not,
    // stop now.
    check_time(noise_batch);
  }

  //============================================================================
  // Cleanup

  // Stop timer
  simulation_timer.stop();

  // Write convergence time
  out.write("\n Time spent converging         : " +
            std::to_string(convergence_timer.elapsed_time()) + " seconds.");
  // Write power iteration time
  out.write("\n Time spent in power iteration : " +
            std::to_string(power_iteration_timer.elapsed_time()) + " seconds.");
  // Write noise time
  out.write("\n Time spent in noise transport : " +
            std::to_string(noise_timer.elapsed_time()) + " seconds.");
  // Write cancellation time
  out.write("\n Time spent in cancellation    : " +
            std::to_string(cancellation_timer.elapsed_time()) + " seconds.\n");

  // Total simulation time
  out.write("\n Total Simulation Time: " +
            std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  // Write flux file
  if (settings::converged && tallies->generations() > 0) {
    tallies->write_tallies();
  }
}

void Noise::power_iteration(bool sample_noise) {
  std::vector<BankedParticle> next_gen;

  if (sample_noise == false) {
    next_gen = transporter->transport(bank, false, nullptr, nullptr);
  } else {
    next_gen = transporter->transport(bank, false, &noise_bank, &noise_maker);
  }
  std::sort(next_gen.begin(), next_gen.end());
  mpi::synchronize();

  if (next_gen.size() == 0) {
    fatal_error("No fission neutrons were produced.");
  }

  // Synchronize the fission banks across all nodes, so that each node will
  // now have the particles that it will be responsible for.
  sync_banks(mpi::node_nparticles, next_gen);

  // Get new keff
  tallies->calc_gen_values();

  // Zero tallies for next generation
  tallies->clear_generation();

  // Do all Pre-Cancelation entropy calculations
  compute_pre_cancellation_entropy(next_gen);

  // Do weight cancelation
  if (settings::regional_cancellation) {
    perform_regional_cancellation(mpi::node_nparticles, next_gen);
  }

  // Calculate net positive and negative weight
  normalize_weights(next_gen);

  // Do all Post-Cancelation entropy calculations
  compute_post_cancellation_entropy(next_gen);

  // Clear and switch particle banks
  bank.clear();

  // Make sure we have the proper history_counter values at each node
  histories_counter = global_histories_counter;
  for (int lower = 0; lower < mpi::rank; lower++) {
    histories_counter += mpi::node_nparticles[static_cast<std::size_t>(lower)];
  }

  bank.reserve(next_gen.size());
  for (auto& p : next_gen) {
    bank.push_back(Particle(p.r, p.u, p.E, p.wgt, histories_counter++));
    bank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
  }

  global_histories_counter += static_cast<uint64_t>(std::accumulate(
      mpi::node_nparticles.begin(), mpi::node_nparticles.end(), 0));

  next_gen.clear();
  next_gen.shrink_to_fit();

  // Output
  pi_generation_output();

  // Zero all entropy bins
  zero_entropy();

  // If we sampled a noise source, sort that too !
  if (sample_noise) std::sort(noise_bank.begin(), noise_bank.end());
}

inline void Noise::perform_regional_cancellation(
    std::vector<uint64_t>& nums, std::vector<BankedParticle>& bank) {
  // To perform cancellation with MPI, what we first do is send ALL of the
  // banked particles back to master. Yes, I know, this is very inefficient. But
  // cancellation works best with the highest density of particles possbile.
  // The smpilest option to make this work without domain decomposition is
  // sending all banked particles to master. Sending all particles to master
  // also avoids some troubles with random number generator seeds, as we can
  // just use the global settings::rng for performing cancellation.
  mpi::Gatherv(bank, 0);

  if (mpi::rank != 0) bank.clear();

  if (mpi::rank == 0) {
    // Only perform cancellation if we are master !!
    std::size_t n_lost_boys = 0;

    // Distribute Particles to Cancellation Bins
    for (auto& p : bank) {
      if (!cancelator->add_particle(p)) n_lost_boys++;
    }

    if (n_lost_boys > 0)
      std::cout << " There are " << n_lost_boys
                << " particles with no cancellation bin.\n";

    // Perform Cancellation for each Bin
    cancelator->perform_cancellation(settings::rng);

    // All particles which were placed into a cancellation bin from next_gen
    // now have modified weights.
    // Now we can get the uniform particles
    auto tmp = cancelator->get_new_particles(settings::rng);
    bank.insert(bank.end(), tmp.begin(), tmp.end());
    nums[0] += tmp.size();

    // All done ! Clear cancelator for next run
    cancelator->clear();
  }

  // Since we added some "uniform particles" to the global fission bank, we need
  // to re-calculate how many particles each node should have, and re-sync all
  // of the banks
  sync_banks(nums, bank);
}

void Noise::noise_simulation() {
  // First, we need to synchronize the initial noise souce bank
  sync_banks(mpi::node_nparticles_noise, noise_bank);

  uint64_t N_noise_tot = static_cast<uint64_t>(std::accumulate(
      mpi::node_nparticles_noise.begin(), mpi::node_nparticles_noise.end(), 0));

  // Start ticking noise timer
  noise_timer.start();

  Output& out = Output::instance();
  out.write("\n -----------------------------------------------\n");

  // Need to keep copy of original so we can override garbage copy which
  // will be produced in this noise batch.
  double original_kcol = tallies->kcol();

  double avg_wgt_mag = 1.;
  if (settings::normalize_noise_source) {
    // Normalize particle weights by the average weight magnitude.
    double sum_wgt_mag = 0.;
    for (const auto& p : noise_bank) {
      sum_wgt_mag += std::sqrt(p.wgt * p.wgt + p.wgt2 * p.wgt2);
    }

    mpi::Allreduce_sum(sum_wgt_mag);

    avg_wgt_mag = sum_wgt_mag / static_cast<double>(N_noise_tot);
    for (auto& p : noise_bank) {
      p.wgt /= avg_wgt_mag;
      p.wgt2 /= avg_wgt_mag;
    }
  }

  // We score the noise source here, because have have normalized
  // all the weights by avg_wgt_mag, so that way when we call
  // rectord_generation on tallies with a multiplier of avg_wgt_mag,
  // we get the correct answer (as if we hadn't normalized noise
  // particle weights).
  tallies->score_noise_source(noise_bank, settings::converged);

  // Cancellation may be performed on noise_bank here

  // Make sure we have the proper history_counter values at each node
  histories_counter = global_histories_counter;
  for (int lower = 0; lower < mpi::rank; lower++) {
    histories_counter +=
        mpi::node_nparticles_noise[static_cast<std::size_t>(lower)];
  }

  // Bank for noise particles.
  std::vector<Particle> nbank;
  for (auto& p : noise_bank) {
    nbank.emplace_back(p.r, p.u, p.E, p.wgt, p.wgt2, histories_counter++);
    nbank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
  }
  global_histories_counter += static_cast<uint64_t>(std::accumulate(
      mpi::node_nparticles_noise.begin(), mpi::node_nparticles_noise.end(), 0));
  noise_bank.clear();

  int noise_gen = 0;
  while (N_noise_tot != 0) {
    noise_gen += 1;

    std::stringstream output;
    output << "  -- " << noise_gen << " running " << N_noise_tot
           << " particles \n";

    auto fission_bank = transporter->transport(nbank, true, nullptr, nullptr);
    if (noise_gen == 1) transported_histories += bank.size();
    std::sort(fission_bank.begin(), fission_bank.end());
    mpi::synchronize();

    // Synchronize the fission banks across all nodes, so that each node will
    // now have the particles that it will be responsible for.
    sync_banks(mpi::node_nparticles_noise, fission_bank);

    // Cancellation may be performed on fission_bank here if we have provided
    // a Cancelator instance.
    if (settings::regional_cancellation_noise &&
        noise_gen <= settings::n_cancel_noise_gens) {
      noise_timer.stop();
      cancellation_timer.start();
      perform_regional_cancellation(mpi::node_nparticles_noise, fission_bank);
      cancellation_timer.stop();
      noise_timer.start();
    }

    nbank.clear();

    // Make sure we have the proper history_counter values at each node
    histories_counter = global_histories_counter;
    for (int lower = 0; lower < mpi::rank; lower++) {
      histories_counter +=
          mpi::node_nparticles_noise[static_cast<std::size_t>(lower)];
    }

    nbank.reserve(fission_bank.size());
    for (const auto& p : fission_bank) {
      nbank.emplace_back(p.r, p.u, p.E, p.wgt, p.wgt2, histories_counter++);
      nbank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
    }
    global_histories_counter += static_cast<uint64_t>(
        std::accumulate(mpi::node_nparticles_noise.begin(),
                        mpi::node_nparticles_noise.end(), 0));

    out.write(output.str());

    // Get the new N_noise_tot
    N_noise_tot = static_cast<std::uint64_t>(
        std::accumulate(mpi::node_nparticles_noise.begin(),
                        mpi::node_nparticles_noise.end(), 0));
  }

  // Get new values
  tallies->calc_gen_values();

  // Keep values
  tallies->record_generation(avg_wgt_mag);

  // Zero tallies for next generation
  tallies->clear_generation();

  // Clear particle banks
  nbank.clear();

  // Output
  noise_output();

  // So that a garbage kcol value isn't sent into PI and screws up the number
  // of particles to make a fissions.
  tallies->set_kcol(original_kcol);

  out.write(" -----------------------------------------------\n\n");

  noise_timer.stop();
}

void Noise::pi_generation_output() {
  std::stringstream output;

  double kcol = tallies->kcol();

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(settings::ngenerations).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << pi_gen << "   " << std::fixed << std::setprecision(5);
  output << kcol;
  output << "   ";

  if (!settings::regional_cancellation && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4);
    output << t_pre_entropy->calculate_entropy();
  }

  if (settings::regional_cancellation && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4)
           << p_pre_entropy->calculate_entropy() << "   ";
    output << n_pre_entropy->calculate_entropy() << "   ";
    output << t_pre_entropy->calculate_entropy() << "   ";

    output << p_post_entropy->calculate_entropy() << "   ";
    output << n_post_entropy->calculate_entropy() << "   ";
    output << t_post_entropy->calculate_entropy() << "   ";
  }

  if (settings::regional_cancellation) {
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Nnet << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Ntot << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Npos << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Nneg << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Wnet << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Wtot << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Wpos << "   ";
    output << std::setw(7) << std::setfill(' ') << std::right;
    output << Wneg;
  }

  // Add line underneath
  output << "\n";

  Output::instance().write(output.str());
}

void Noise::noise_output() {
  std::stringstream output;

  output << "\n";

  // Write the current noise batch
  output << " Simulated noise batch " << noise_batch << "/"
         << settings::ngenerations << ".";

  output << "\n";

  Output::instance().write(output.str());
}

void Noise::normalize_weights(std::vector<BankedParticle>& next_gen) {
  double W = 0.;
  double W_neg = 0.;
  double W_pos = 0.;
  Nnet = 0;
  Ntot = 0;
  Npos = 0;
  Nneg = 0;
  for (size_t i = 0; i < next_gen.size(); i++) {
    if (next_gen[i].wgt > 0.) {
      W_pos += next_gen[i].wgt;
      Npos++;
    } else {
      W_neg -= next_gen[i].wgt;
      Nneg++;
    }
  }

  mpi::Allreduce_sum(W);
  mpi::Allreduce_sum(W_neg);
  mpi::Allreduce_sum(W_pos);
  mpi::Allreduce_sum(Npos);
  mpi::Allreduce_sum(Nneg);

  W = W_pos - W_neg;
  Ntot = Npos + Nneg;
  Nnet = Npos - Nneg;

  // Re-Normalize particle weights
  double w_per_part = static_cast<double>(settings::nparticles) / W;
  W *= w_per_part;
  W_neg *= w_per_part;
  W_pos *= w_per_part;

  for (std::size_t i = 0; i < next_gen.size(); i++) {
    next_gen[i].wgt *= w_per_part;
  }

  double Wtt = W_pos + W_neg;
  Wnet = static_cast<int>(std::round(W));
  Wtot = static_cast<int>(std::round(Wtt));
  Wpos = static_cast<int>(std::round(W_pos));
  Wneg = static_cast<int>(std::round(W_neg));
}

void Noise::premature_kill() {
  // See if the user really wants to kill the program
  bool user_said_kill = false;

  if (mpi::rank == 0 && terminate == false) {
    std::string response;
    std::cout << "\n Do you really want to stop the simulation ? (y/N) => ";
    std::cin >> response;
    if (response == "y" || response == "Y") user_said_kill = true;
  }

  mpi::Bcast<bool>(user_said_kill, 0);

  if (user_said_kill || terminate) {
    Output::instance().write(
        "\n Simulation has been stopped by user. Cleaning up...");
    // Looks like they really want to kill it.

    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = noise_batch;
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

void Noise::compute_pre_cancellation_entropy(
    std::vector<BankedParticle>& next_gen) {
  if (t_pre_entropy && settings::regional_cancellation) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      p_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      n_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      t_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }
    p_pre_entropy->synchronize_entropy_across_nodes();
    n_pre_entropy->synchronize_entropy_across_nodes();
    t_pre_entropy->synchronize_entropy_across_nodes();
  } else if (t_pre_entropy) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      t_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }
    t_pre_entropy->synchronize_entropy_across_nodes();
  }
}

void Noise::compute_post_cancellation_entropy(
    std::vector<BankedParticle>& next_gen) {
  if (t_pre_entropy && settings::regional_cancellation) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      p_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      n_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      t_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }
    p_post_entropy->synchronize_entropy_across_nodes();
    n_post_entropy->synchronize_entropy_across_nodes();
    t_post_entropy->synchronize_entropy_across_nodes();
  }
}

void Noise::zero_entropy() {
  if (p_pre_entropy) {
    p_pre_entropy->zero();
  }

  if (n_pre_entropy) {
    n_pre_entropy->zero();
  }

  if (t_pre_entropy) {
    t_pre_entropy->zero();
  }

  if (p_post_entropy) {
    p_post_entropy->zero();
  }

  if (n_post_entropy) {
    n_post_entropy->zero();
  }

  if (t_post_entropy) {
    t_post_entropy->zero();
  }
}
