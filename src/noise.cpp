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
#include <utils/kahan.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <xtensor/xtensor.hpp>

#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

// Here, we sample a particle source exactly as we do for power iteration.
void Noise::initialize() {
  this->write_output_info();

  if (in_source_file_name.empty() && sources.empty()) {
    fatal_error("No sources or source file provided.");
  }

  // If not using source file
  if (in_source_file_name.empty()) {
    sample_source_from_sources();
  } else {
    // Load source file
    load_source_from_file();
  }

  mpi::node_nparticles_noise.resize(static_cast<std::size_t>(mpi::size), 0);

  Tallies::instance().allocate_batch_arrays(nbatches * nskip + nignored);
  Tallies::instance().verify_track_length_tallies(
      noise_particle_mover->track_length_compatible());
}

void Noise::write_output_info() const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  h5.createAttribute<std::string>("simulation-mode", "noise");
  h5.createAttribute("nparticles", nparticles);
  h5.createAttribute("nbatches", nbatches);
  h5.createAttribute("nignored", nignored);
  h5.createAttribute("nskip", nskip);
  h5.createAttribute("combing", combing);
  h5.createAttribute("ncancel-noise-gens", ncancel_noise_gens);
  h5.createAttribute("normalize-noise-source", normalize_noise_source_);
  h5.createAttribute("regional-cancellation", regional_cancellation_);
  h5.createAttribute("noise-regional-cancellation",
                     regional_cancellation_noise_);
  if (in_source_file_name.empty() == false)
    h5.createAttribute("source-file", in_source_file_name);

  if (regional_cancellation_ || regional_cancellation_noise_) {
    auto cancelator_grp = h5.createGroup("cancelator");
    cancelator->write_output_info(cancelator_grp);
  }

  noise_particle_mover->write_output_info("noise-");
  particle_mover->write_output_info("power-iterator-");
}

void Noise::set_entropy(const Entropy& entropy) {
  t_pre_entropy = std::make_shared<Entropy>(entropy);
  t_pre_entropy->set_sign(Entropy::Sign::Total);

  if (cancelator) {
    p_pre_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    p_pre_entropy->set_sign(Entropy::Sign::Positive);

    n_pre_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    n_pre_entropy->set_sign(Entropy::Sign::Negative);

    t_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    t_post_entropy->set_sign(Entropy::Sign::Total);

    p_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    p_post_entropy->set_sign(Entropy::Sign::Positive);

    n_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    n_post_entropy->set_sign(Entropy::Sign::Negative);
  }
}

void Noise::set_cancelator(std::shared_ptr<Cancelator> cncl) {
  cancelator = cncl;
  if (cancelator == nullptr) {
    fatal_error("Was given a nullptr for a Cancelator.");
  }

  // If we already had an entropy, we need to copy it to all others now
  if (t_pre_entropy) {
    p_pre_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    p_pre_entropy->set_sign(Entropy::Sign::Positive);

    n_pre_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    n_pre_entropy->set_sign(Entropy::Sign::Negative);

    t_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    t_post_entropy->set_sign(Entropy::Sign::Total);

    p_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    p_post_entropy->set_sign(Entropy::Sign::Positive);

    n_post_entropy = std::make_shared<Entropy>(*t_pre_entropy);
    n_post_entropy->set_sign(Entropy::Sign::Negative);
  }
}

void Noise::add_source(std::shared_ptr<Source> src) {
  if (src)
    sources.push_back(src);
  else
    fatal_error("Was provided a nullptr Source.");
}

void Noise::set_regional_cancellation(bool rc) {
  if (rc && cancelator == nullptr) {
    fatal_error("Cannot use regional cancellation without a Cancelator.");
  }

  regional_cancellation_ = rc;
}

void Noise::set_regional_cancellation_noise(bool rcn) {
  if (rcn && cancelator == nullptr) {
    fatal_error(
        "Cannot do regional cancellation on noise without a Cancelator.");
  }

  regional_cancellation_noise_ = rcn;
}

void Noise::load_source_from_file() {
  // Load source file
  auto h5s = H5::File(in_source_file_name, H5::File::ReadOnly);

  // Get the dataset and dimensions
  auto source_ds = h5s.getDataSet("source");
  std::vector<std::size_t> dimensions = source_ds.getSpace().getDimensions();

  // Check dimensions
  if (dimensions.size() != 2 || dimensions[1] != 9) {
    fatal_error("Invalid source from file dimensions.");
  }

  // Read in array
  xt::xtensor<double, 2> source;
  source.resize(dimensions);
  source_ds.read<double>(source.data());

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

  // Distribute the remainder particles amongst the nodes. There are at most
  // mpi::size-1 remainder particles, so we will distribute them until we
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
  nparticles = static_cast<std::size_t>(std::round(tot_wgt));
  Tallies::instance().set_total_weight(std::round(tot_wgt));
}

void Noise::sample_source_from_sources() {
  Output::instance().write(" Generating source particles...\n");
  // Calculate the base number of particles per node to run
  uint64_t base_particles_per_node =
      static_cast<uint64_t>(nparticles / static_cast<std::size_t>(mpi::size));
  uint64_t remainder =
      static_cast<uint64_t>(nparticles % static_cast<std::size_t>(mpi::size));

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
  bank = sample_sources(
      sources, mpi::node_nparticles[static_cast<std::size_t>(mpi::rank)]);

  global_histories_counter += static_cast<uint64_t>(std::accumulate(
      mpi::node_nparticles.begin(), mpi::node_nparticles.end(), 0));

  Tallies::instance().set_total_weight(static_cast<double>(nparticles));
}

void Noise::print_header() const {
  std::stringstream output;
  output << "\n";

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(nbatches * nskip + nignored).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << "Gen"
         << "      k     ";

  if (regional_cancellation_ == false && t_pre_entropy) {
    output << "    Entropy  ";
  }

  // Add line underneath
  output << "\n " << std::string(output.str().size() - 3, '-') << "\n";

  Output::instance().write(output.str());
}

bool Noise::out_of_time(std::size_t batch) {
  // Get the average time per batch
  double T_avg = noise_batch_timer.elapsed_time() / static_cast<double>(batch);

  // See how much time we have used so far.
  double T_used = settings::alpha_omega_timer.elapsed_time();

  // See how much time we have left
  double T_remaining = settings::max_time - T_used;

  // If the remaining time is less than 3*T_avg, than we just kill it
  // it here, so that we are sure we don't run over.
  if (T_remaining < 3. * T_avg) {
    return true;
  }

  return false;
}

void Noise::check_time(std::size_t batch) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(batch);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting nbatches to batch will cause us to exit the
    // transport loop as if we finished the simulation normally.
    nbatches = batch;
  }
}

void Noise::run() {
  Output& out = Output::instance();
  out.write(" Running Noise Problem...\n");
  out.write(" NPARTICLES: " + std::to_string(nparticles) + ", ");
  out.write(" NBATCHES: " + std::to_string(nbatches) + ",");
  out.write(" NIGNORED: " + std::to_string(nignored) + ",");
  out.write(" NSKIP: " + std::to_string(nskip) + "\n\n");
  print_header();

  // Zero all entropy bins
  zero_entropy();

  // Start timer
  simulation_timer.reset();
  mpi::synchronize();  // Make sure everyone is here before starting the timer
  simulation_timer.start();

  //============================================================================
  // Run standard power iteration until the source has converged.
  // Make sure converged is false so we don't collect statistics.
  convergence_timer.reset();
  convergence_timer.start();
  converged = false;
  Tallies::instance().set_scoring(converged);
  for (pi_gen = 1; pi_gen <= nignored; pi_gen++) {
    power_iteration(false);
  }

  // We now consider the source to be converged.
  converged = true;
  Tallies::instance().set_scoring(converged);

  // Make sure timers is cleared
  convergence_timer.stop();
  noise_timer.reset();
  cancellation_timer.reset();
  power_iteration_timer.reset();

  //============================================================================
  // Now we run the desired number of noise batches
  noise_batch_timer.reset();
  noise_batch_timer.start();
  for (noise_batch = 1; noise_batch <= nbatches; noise_batch++) {
    // We first run the decorrelation PI batches, where no source is sampled.
    // This should be nskip - 1.
    // Make sure converged is false so we don't collect statistics on the
    // power iteration parts.
    power_iteration_timer.start();
    for (int skip = 0; skip < static_cast<int>(nskip) - 1; skip++) {
      power_iteration(false);
      pi_gen++;
    }

    // We now run the last remaining power iteration, but here we sample the
    // noise source
    power_iteration(true);
    power_iteration_timer.stop();
    pi_gen++;

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
  if (converged && Tallies::instance().generations() > 0) {
    Tallies::instance().write_tallies(
        particle_mover->track_length_compatible());
  }
}

void Noise::power_iteration(bool sample_noise) {
  // We set converged to false so that we don't waste time doing tallies in the
  // power iteration portion of the simulation. We save a copy to reset after.
  const bool orig_converged = converged;
  converged = false;
  Tallies::instance().set_scoring(converged);

  std::vector<BankedParticle> next_gen;
  if (sample_noise == false) {
    next_gen = particle_mover->transport(bank, std::nullopt, nullptr, nullptr);
  } else {
    next_gen = particle_mover->transport(bank, std::nullopt, &noise_bank,
                                         &noise_maker);
  }

  if (next_gen.size() == 0) {
    fatal_error("No fission neutrons were produced.");
  }

  // Do all Pre-Cancelation entropy calculations
  compute_pre_cancellation_entropy(next_gen);

  // Get new keff
  Tallies::instance().calc_gen_values();

  // Zero tallies for next generation
  Tallies::instance().clear_generation();

  // Do weight cancelation
  if (regional_cancellation_) {
    perform_regional_cancellation(cancelator, next_gen);
  }

  // Calculate net positive and negative weight
  normalize_weights(nparticles, next_gen, Npos, Nneg, Nnet, Ntot, Wpos, Wneg,
                    Wnet, Wtot);

  // Comb particles
  if (combing) comb_particles(next_gen, Npos, Nneg, Nnet, Ntot);

  // Synchronize the fission banks across all nodes, so that each node will
  // now have the particles that it will be responsible for.
  sync_banks(mpi::node_nparticles, next_gen);

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

  // Reset converged to original value
  converged = orig_converged;
  Tallies::instance().set_scoring(converged);
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
  double original_kcol = Tallies::instance().kcol();

  double avg_wgt_mag = 1.;
  if (normalize_noise_source_) {
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
  Tallies::instance().score_noise_source(noise_bank);

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

  std::size_t noise_gen = 0;
  while (N_noise_tot != 0) {
    noise_gen += 1;

    std::stringstream output;
    output << "  -- " << noise_gen << " running " << N_noise_tot
           << " particles \n";

    auto fission_bank =
        noise_particle_mover->transport(nbank, noise_params, nullptr, nullptr);
    if (noise_gen == 1) transported_histories += bank.size();

    // Cancellation may be performed on fission_bank here if we have provided
    // a Cancelator instance.
    if (regional_cancellation_noise_ && noise_gen <= ncancel_noise_gens) {
      noise_timer.stop();
      cancelator->set_cancel_dual_weights(true);
      cancellation_timer.start();
      perform_regional_cancellation(cancelator, fission_bank);
      cancellation_timer.stop();
      cancelator->set_cancel_dual_weights(false);
      noise_timer.start();
    }

    // Synchronize the fission banks across all nodes, so that each node will
    // now have the particles for which it will be responsible.
    sync_banks(mpi::node_nparticles_noise, fission_bank);

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
  Tallies::instance().calc_gen_values();

  // Keep values
  Tallies::instance().record_generation(avg_wgt_mag);

  // Zero tallies for next generation
  Tallies::instance().clear_generation();

  // Clear particle banks
  nbank.clear();

  // Output
  noise_output();

  // So that a garbage kcol value isn't sent into PI and screws up the number
  // of particles to make a fissions.
  Tallies::instance().set_kcol(original_kcol);

  out.write(" -----------------------------------------------\n\n");

  noise_timer.stop();
}

void Noise::pi_generation_output() {
  std::stringstream output;

  double kcol = Tallies::instance().kcol();

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(nbatches * nskip + nignored).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << pi_gen << "   " << std::fixed << std::setprecision(5);
  output << kcol;
  output << "   ";

  if (regional_cancellation_ == false && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4);
    output << t_pre_entropy->calculate_entropy();
  }

  if (regional_cancellation_ && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4)
           << p_pre_entropy->calculate_entropy() << "   ";
    output << n_pre_entropy->calculate_entropy() << "   ";
    output << t_pre_entropy->calculate_entropy() << "   ";

    output << p_post_entropy->calculate_entropy() << "   ";
    output << n_post_entropy->calculate_entropy() << "   ";
    output << t_post_entropy->calculate_entropy() << "   ";
  }

  if (regional_cancellation_) {
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
  output << " Simulated noise batch " << noise_batch << "/" << nbatches << ".";

  output << "\n";

  Output::instance().write(output.str());
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
    nbatches = noise_batch;
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

void Noise::compute_pre_cancellation_entropy(
    std::vector<BankedParticle>& next_gen) {
  if (t_pre_entropy && regional_cancellation_) {
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
  if (t_pre_entropy && regional_cancellation_) {
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

std::shared_ptr<Noise> make_noise_simulator(const YAML::Node& sim) {
  // Get the number of particles
  if (!sim["nparticles"] || sim["nparticles"].IsScalar() == false) {
    fatal_error("No nparticles entry in noise simulation.");
  }
  std::size_t nparticles = sim["nparticles"].as<std::size_t>();

  // Get the number of noise batches
  if (!sim["nbatches"] || sim["nbatches"].IsScalar() == false) {
    fatal_error("No nbatches entry in noise simulation.");
  }
  std::size_t ngenerations = sim["nbatches"].as<std::size_t>();

  // Get the number of ignored generations
  if (!sim["nignored"] || sim["nignored"].IsScalar() == false) {
    fatal_error("No nignored entry in noise simulation.");
  }
  std::size_t nignored = sim["nignored"].as<std::size_t>();

  // Get the number of skipped generations between noise sampling
  if (!sim["nskip"] || sim["nskip"].IsScalar() == false) {
    fatal_error("No nskip entry in noise simulation.");
  }
  std::size_t nskip = sim["nskip"].as<std::size_t>();

  std::string in_source_file;
  if (sim["source-file"] && sim["source-file"].IsScalar()) {
    in_source_file = sim["source-file"].as<std::string>();
  } else if (sim["source-file"]) {
    fatal_error("Invalid source-file entry in k-eigenvalue simulation.");
  }

  // Read all sources if we don't have a source file
  std::vector<std::shared_ptr<Source>> sources;
  if (in_source_file.empty()) {
    if (!sim["sources"] || sim["sources"].IsSequence() == false) {
      fatal_error("No sources entry in noise simulation.");
    }
    for (std::size_t s = 0; s < sim["sources"].size(); s++) {
      sources.push_back(make_source(sim["sources"][s]));
    }
  }

  // Combing
  bool combing = false;
  if (sim["combing"] && sim["combing"].IsScalar()) {
    combing = sim["combing"].as<bool>();
  } else if (sim["combing"]) {
    fatal_error("Invalid combing entry in noise simulation.");
  }

  // For noise, we need a NoiseMaker. We will now read in all of the provided
  // noise sources in the simulation entry
  NoiseMaker noise_maker;
  if (!sim["noise-sources"] || sim["noise-sources"].IsSequence() == false) {
    fatal_error("No noise-sources entry in noise simulation.");
  }
  for (std::size_t s = 0; s < sim["noise-sources"].size(); s++) {
    noise_maker.add_noise_source(sim["noise-sources"][s]);
  }

  bool normalize_noise_source = true;
  if (sim["normalize-noise-source"] &&
      sim["normalize-noise-source"].IsScalar()) {
    normalize_noise_source = sim["normalize-noise-source"].as<bool>();
  } else if (sim["normalize-noise-source"]) {
    fatal_error(
        "Invalid normalize-noise-source entry provided in noise simulation.");
  }

  // Get entropy bins
  std::optional<Entropy> entropy = std::nullopt;
  if (sim["entropy"] && sim["entropy"].IsMap()) {
    entropy = make_entropy(sim["entropy"]);
  } else if (sim["entropy"]) {
    fatal_error("Invalid entropy entry in noise simulation.");
  }

  // For noise simulations, we need two different particle movers. One for the
  // noise particles and one for the normal particles in power iteration. We
  // start by making the particle mover for noise. It will use the prescribed
  // transport operator, but will always use NoiseBranchingCollision.
  std::shared_ptr<INoiseParticleMover> noise_mover =
      make_noise_particle_mover(sim);

  // Now we can get the particle mover for power iteration. This will sue the
  // same transport operator, but will use the prescribed collision operator
  // in the input file.
  std::shared_ptr<IParticleMover> pi_mover =
      make_particle_mover(sim, settings::SimMode::KEFF);

  std::size_t ncancel_noise_gens = std::numeric_limits<std::size_t>::max();
  if (sim["ncancel-noise-gens"] && sim["ncancel-noise-gens"].IsScalar()) {
    ncancel_noise_gens = sim["ncancel-noise-gens"].as<std::size_t>();
  } else if (sim["ncancel-noise-gens"]) {
    fatal_error(
        "Invalid ncancel-noise-gens entry provided in noise simulation.");
  }

  // Get cancelator
  std::shared_ptr<Cancelator> cancelator = nullptr;
  if (sim["cancelator"] && sim["cancelator"].IsMap()) {
    cancelator = make_cancelator(sim["cancelator"]);
  } else if (sim["cancelator"]) {
    fatal_error("Invalid cancelator entry in noise simulation.");
  }

  // Check for regional cancellation options. By default, we assume noise
  // cancellation if a cancelator is present.
  bool noise_cancellation = cancelator ? true : false;
  if (sim["noise-cancellation"] && sim["noise-cancellation"].IsScalar()) {
    noise_cancellation = sim["noise-cancellation"].as<bool>();
  } else if (sim["noise-cancellation"]) {
    fatal_error("Invalid noise-cancellation entry in noise simulation.");
  }

  bool cancellation = false;
  if (sim["cancellation"] && sim["cancellation"].IsScalar()) {
    noise_cancellation = sim["cancellation"].as<bool>();
  } else if (sim["cancellation"]) {
    fatal_error("Invalid cancellation entry in noise simulation.");
  }

  if ((noise_cancellation || cancellation) && cancelator == nullptr) {
    fatal_error(
        "Cancellation is desired, but no cancelator was provided in noise "
        "simulation.");
  }

  // Make sure our transport operators are compatible with cancellation !
  if (noise_cancellation) {
    cancelator->check_particle_mover_compatibility(noise_mover);
  }
  if (cancellation) {
    cancelator->check_particle_mover_compatibility(pi_mover);
  }

  // Read in noise parameters
  double omega;
  if (!sim["noise-angular-frequency"] ||
      sim["noise-angular-frequency"].IsScalar() == false) {
    fatal_error(
        "No valid noise-angular-frequency entry provided in noise simulation.");
  }
  omega = sim["noise-angular-frequency"].as<double>();

  double keff;
  if (!sim["keff"] || sim["keff"].IsScalar() == false) {
    fatal_error("No valid keff entry provided in noise simulation.");
  }
  keff = sim["keff"].as<double>();

  double eta = 1.;
  if (sim["eta"] && sim["eta"].IsScalar()) {
    eta = sim["eta"].as<double>();
  } else if (sim["eta"]) {
    fatal_error("Invalid eta entry provided in noise simulation.");
  }

  NoiseParameters noise_params;
  noise_params.omega = omega;
  noise_params.keff = keff;
  noise_params.eta = eta;

  noise_maker.noise_parameters() = noise_params;

  // Create simulation
  std::shared_ptr<Noise> simptr =
      std::make_shared<Noise>(noise_mover, pi_mover, noise_params, noise_maker);

  // Set quantities
  simptr->set_nparticles(nparticles);
  simptr->set_nbatches(ngenerations);
  simptr->set_nignored(nignored);
  simptr->set_nskip(nskip);
  simptr->set_combing(combing);
  simptr->set_normalize_noise_source(normalize_noise_source);
  simptr->set_ncancel_noise_gens(ncancel_noise_gens);
  if (cancelator) simptr->set_cancelator(cancelator);
  simptr->set_regional_cancellation(cancellation);
  simptr->set_regional_cancellation_noise(noise_cancellation);
  for (const auto& src : sources) {
    simptr->add_source(src);
  }
  if (in_source_file.empty() == false) {
    simptr->set_in_source_file(in_source_file);
  }
  if (entropy) simptr->set_entropy(*entropy);

  return simptr;
}
