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
#include <geometry/geometry.hpp>
#include <simulation/power_iterator.hpp>
#include <utils/error.hpp>
#include <utils/kahan.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <xtensor/xtensor.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <vector>

void PowerIterator::initialize() {
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

  // Assign History ID as family ID. At beginning, history IDs should start at
  // 1, so this SHOULD yield the desired results.
  for (auto& p : bank) p.set_family_id(p.history_id());

  Tallies::instance().allocate_batch_arrays(ngenerations);
  Tallies::instance().verify_track_length_tallies(
      particle_mover->track_length_compatible());
}

void PowerIterator::write_output_info() const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  h5.createAttribute<std::string>("simulation-mode", "k-eigenvalue");
  h5.createAttribute("nparticles", nparticles);
  h5.createAttribute("ngenerations", ngenerations);
  h5.createAttribute("nignored", nignored);
  h5.createAttribute("combing", combing);
  h5.createAttribute("families", calc_families);
  h5.createAttribute("pair-distance-sqrd", pair_distance_sqrd);
  h5.createAttribute("empty-entropy-bins-frac", empty_entropy_bins);

  if (in_source_file_name.empty() == false)
    h5.createAttribute("source-file", in_source_file_name);

  if (cancelator) {
    h5.createAttribute("regional-cancellation", true);
    auto cancelator_grp = h5.createGroup("cancelator");
    cancelator->write_output_info(cancelator_grp);
  } else {
    h5.createAttribute("regional-cancellation", false);
  }

  if (t_pre_entropy) {
    h5.createAttribute("entropy", true);
  } else {
    h5.createAttribute("entropy", true);
  }

  particle_mover->write_output_info();
}

void PowerIterator::set_cancelator(std::shared_ptr<Cancelator> cncl) {
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

void PowerIterator::set_entropy(const Entropy& entropy) {
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

void PowerIterator::add_source(std::shared_ptr<Source> src) {
  if (src)
    sources.push_back(src);
  else
    fatal_error("Was provided a nullptr Source.");
}

void PowerIterator::load_source_from_file() {
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
        mpi::node_nparticles[static_cast<std::size_t>(lower_rank)];
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

void PowerIterator::sample_source_from_sources() {
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

void PowerIterator::print_header() {
  Output& out = Output::instance();
  out.write("\n");
  std::stringstream output;

  // First get the number of columns required to print the max generation number
  int n_col_gen = static_cast<int>(std::to_string(ngenerations).size());

  if (cancelator && t_pre_entropy) {
    out.write(
        std::string(static_cast<std::size_t>(std::max(n_col_gen, 3) + 1), ' ') +
        "                                         Pre-Cancelation Entropies    "
        "         Post-Cancelation Entropies\n");
    out.write(
        std::string(static_cast<std::size_t>(std::max(n_col_gen, 3) + 1), ' ') +
        "                                   "
        "------------------------------------   "
        "------------------------------------\n");
  }

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << "Gen"
         << "      k         kavg +/- err    ";

  if (cancelator == nullptr && t_pre_entropy) {
    output << "    Entropy  ";
  }

  if (cancelator && t_pre_entropy) {
    output << "     Positive     Nevative        Total     Positive     "
              "Negative        Total";
  }

  if (cancelator) {
    output << "      Nnet      Ntot      Npos      Nneg      Wnet      Wtot    "
              "  Wpos      Wneg";
  }

  if (pair_distance_sqrd) {
    output << "      <r^2>  ";
  }

  if (calc_families) {
    output << "   N_Families  ";
  }

  if (t_pre_entropy && empty_entropy_bins) {
    output << "   Empty Entropy Frac.  ";
  }

  // Add line underneath
  output << "\n " << std::string(output.str().size() - 3, '-') << "\n";

  out.write(output.str());
}

void PowerIterator::generation_output() {
  std::stringstream output;

  double kcol = Tallies::instance().kcol();

  if (gen == 1) print_header();

  // First get the number of columns required to print the max generation number
  int n_col_gen = static_cast<int>(std::to_string(ngenerations).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << gen << "   " << std::fixed << std::setprecision(5);
  output << kcol;
  if (gen <= nignored + 1) {
    output << "                      ";
  } else {
    double kcol_avg = Tallies::instance().kcol_avg();
    double kcol_err = Tallies::instance().kcol_err();
    output << "   " << kcol_avg << " +/- " << kcol_err;
  }
  output << "   ";

  if (cancelator == nullptr && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4);
    output << t_pre_entropy_vec.back();
  }

  if (cancelator && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4)
           << p_pre_entropy_vec.back() << "   ";
    output << n_pre_entropy_vec.back() << "   ";
    output << t_pre_entropy_vec.back() << "   ";

    output << p_post_entropy_vec.back() << "   ";
    output << n_post_entropy_vec.back() << "   ";
    output << t_post_entropy_vec.back() << "   ";
  }

  if (cancelator) {
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

  if (pair_distance_sqrd) {
    output << "   " << std::setw(10) << std::scientific << std::setprecision(4)
           << r_sqrd;
  }

  if (calc_families) {
    std::size_t nfamilies = families_set.size();
    mpi::Reduce_sum(nfamilies, 0);
    families_vec.push_back(nfamilies);
    output << "   " << std::fixed << nfamilies;
  }

  if (t_pre_entropy && empty_entropy_bins) {
    output << "   " << std::setprecision(6) << empty_entropy_frac_vec.back();
  }

  // Add line underneath
  output << "\n";

  Output::instance().write(output.str());
}

void PowerIterator::run() {
  Output& out = Output::instance();
  out.write(" Running k Eigenvalue Problem...\n");
  out.write(" NPARTICLES: " + std::to_string(nparticles) + ", ");
  out.write(" NGENERATIONS: " + std::to_string(ngenerations) + ",");
  out.write(" NIGNORED: " + std::to_string(nignored) + "\n");

  // Zero all entropy bins
  zero_entropy();

  // Start timer
  simulation_timer.reset();
  mpi::synchronize();  // Make sure everyone is here before starting the timer
  simulation_timer.start();

  // Make sure we have room for the generation times
  generation_times.reserve(ngenerations);

  // Check for immediate convergence (i.e. we read source from file)
  Tallies::instance().set_scoring(false);
  if (nignored == 0) {
    Tallies::instance().set_scoring(true);
  }

  for (gen = 1; gen <= ngenerations; gen++) {
    generation_timer.start();

    if (calc_families) {
      // Before transport, we get the set of all families, so that we can know
      // how many particle families entered the generation.
      families_set.clear();
      for (auto& p : bank) families_set.insert(p.family_id());
    }

    std::vector<BankedParticle> next_gen = particle_mover->transport(bank);
    transported_histories += bank.size();

    if (next_gen.empty()) {
      fatal_error("No fission neutrons were produced.");
    }

    // Do all Pre-Cancelation entropy calculations
    compute_pre_cancellation_entropy(next_gen);

    // Get new keff
    Tallies::instance().calc_gen_values();

    // Do weight cancelation
    if (cancelator) {
      perform_regional_cancellation(cancelator, next_gen);
    }

    // Calculate net positive and negative weight
    normalize_weights(nparticles, next_gen, Npos, Nneg, Nnet, Ntot, Wpos, Wneg,
                      Wnet, Wtot);
    save_weights();

    // Comb particles
    if (combing) comb_particles(next_gen, Npos, Nneg, Nnet, Ntot);

    // Synchronize particles across nodes
    sync_banks(mpi::node_nparticles, next_gen);

    // Compute pair distance squared
    if (pair_distance_sqrd) compute_pair_dist_sqrd(next_gen);

    // Score the source and gen if passed nignored
    if (gen > nignored) {
      Tallies::instance().score_source(next_gen);
      Tallies::instance().record_generation();
    }

    // Zero tallies for next generation
    Tallies::instance().clear_generation();

    // Do all Post-Cancelation entropy calculations
    compute_post_cancellation_entropy(next_gen);

    // Clear and switch particle banks
    bank.clear();

    // Make sure we have the proper history_counter values at each node
    histories_counter = global_histories_counter;
    for (int lower = 0; lower < mpi::rank; lower++) {
      histories_counter +=
          mpi::node_nparticles[static_cast<std::size_t>(lower)];
    }

    bank.reserve(next_gen.size());
    for (auto& p : next_gen) {
      bank.push_back(Particle(p.r, p.u, p.E, p.wgt, histories_counter++));
      bank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
      bank.back().set_family_id(p.family_id);
    }

    global_histories_counter += static_cast<uint64_t>(std::accumulate(
        mpi::node_nparticles.begin(), mpi::node_nparticles.end(), 0));

    next_gen.clear();
    next_gen.shrink_to_fit();

    // Output
    generation_output();

    // Check if signal has been sent after generation keff has been
    // recorded, and cancellation has occurred. Otherwise source file
    // will be larger than necessary, and wrong number of gens will be
    // in output file based on number averaged for tallies.
    sync_signaled();
    if (signaled) premature_kill();

    // Check if we have enough time to finish the simulation. If not,
    // stop now.
    check_time(gen);

    // Zero all entropy bins
    zero_entropy();

    // Once ignored generations are finished, mark as true to start
    // doing tallies
    if (gen == nignored) {
      Tallies::instance().set_scoring(true);
    }

    generation_timer.stop();
    generation_times.push_back(generation_timer.elapsed_time());
    generation_timer.reset();
  }

  // Stop timer
  simulation_timer.stop();
  out.write("\n Total Simulation Time: " +
            std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  if (gen > nignored) {
    // Write the final results of all estimators
    out.write("\n");
    std::stringstream output;
    output << " Results using " << gen - nignored << " active generations:\n";
    output << " -----------------------------------\n";
    output << std::fixed << std::setprecision(6);
    output << " | kcol    = " << Tallies::instance().kcol_avg() << " +/- "
           << Tallies::instance().kcol_err() << " |\n";
    if (settings::energy_mode == settings::EnergyMode::CE) {
      output << " | kabs    = " << Tallies::instance().kabs_avg() << " +/- "
             << Tallies::instance().kabs_err() << " |\n";
    }
    if (particle_mover->track_length_compatible()) {
      output << " | ktrk    = " << Tallies::instance().ktrk_avg() << " +/- "
             << Tallies::instance().ktrk_err() << " |\n";
    }
    output << " | leakage = " << Tallies::instance().leakage_avg() << " +/- "
           << Tallies::instance().leakage_err() << " |\n";
    output << " -----------------------------------\n";

    out.write(output.str());

    // Write saved warnings
    out.write_saved_warnings();

    // Save other outputs
    out.write("\n");
  }

  // Write results file
  Tallies::instance().write_tallies(particle_mover->track_length_compatible());
  write_entropy_families_etc_to_results();

  // write source
  if (save_source) write_source(bank);
}

void PowerIterator::write_entropy_families_etc_to_results() const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  if (families_vec.size() > 0) {
    h5.createDataSet("results/families", families_vec);
  }

  if (r_sqrd_vec.size() > 0) {
    h5.createDataSet("results/pair-dist-sqrd", r_sqrd_vec);
  }

  if (t_pre_entropy_vec.size() > 0 && cancelator == nullptr) {
    h5.createDataSet("results/entropy", t_pre_entropy_vec);
  } else if (t_pre_entropy_vec.size() > 0) {
    h5.createDataSet("results/total-pre-cancel-entropy", t_pre_entropy_vec);
    h5.createDataSet("results/neg-pre-cancel-entropy", n_pre_entropy_vec);
    h5.createDataSet("results/pos-pre-cancel-entropy", p_pre_entropy_vec);

    h5.createDataSet("results/total-post-cancel-entropy", t_post_entropy_vec);
    h5.createDataSet("results/neg-post-cancel-entropy", n_post_entropy_vec);
    h5.createDataSet("results/pos-post-cancel-entropy", p_post_entropy_vec);
  }

  if (empty_entropy_bins) {
    h5.createDataSet("results/empty-entropy-frac", empty_entropy_frac_vec);
  }

  if (Nnet_vec.size() > 0) {
    h5.createDataSet("results/Nnet", Nnet_vec);
  }

  if (Ntot_vec.size() > 0) {
    h5.createDataSet("results/Ntot", Ntot_vec);
  }

  if (Npos_vec.size() > 0) {
    h5.createDataSet("results/Npos", Npos_vec);
  }

  if (Nneg_vec.size() > 0) {
    h5.createDataSet("results/Nneg", Nneg_vec);
  }

  if (Wnet_vec.size() > 0) {
    h5.createDataSet("results/Wnet", Wnet_vec);
  }

  if (Wtot_vec.size() > 0) {
    h5.createDataSet("results/Wtot", Wtot_vec);
  }

  if (Wpos_vec.size() > 0) {
    h5.createDataSet("results/Wpos", Wpos_vec);
  }

  if (Wneg_vec.size() > 0) {
    h5.createDataSet("results/Wneg", Wneg_vec);
  }

  // Write simulation time and number of particle transported
  h5.createAttribute("simulation-time", simulation_timer.elapsed_time());
  h5.createAttribute("nparticles-transported", histories_counter);
  h5.createAttribute("generation-times", generation_times);
}

void PowerIterator::save_weights() {
  // Save values to vectors
  Nnet_vec.push_back(Nnet);
  Ntot_vec.push_back(Ntot);
  Npos_vec.push_back(Npos);
  Nneg_vec.push_back(Nneg);

  Wnet_vec.push_back(Wnet);
  Wtot_vec.push_back(Wtot);
  Wpos_vec.push_back(Wpos);
  Wneg_vec.push_back(Wneg);
}

void PowerIterator::compute_pre_cancellation_entropy(
    std::vector<BankedParticle>& next_gen) {
  if (t_pre_entropy && cancelator) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      p_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      n_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      t_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }

    p_pre_entropy->synchronize_entropy_across_nodes();
    n_pre_entropy->synchronize_entropy_across_nodes();
    t_pre_entropy->synchronize_entropy_across_nodes();

    p_pre_entropy_vec.push_back(p_pre_entropy->calculate_entropy());
    n_pre_entropy_vec.push_back(n_pre_entropy->calculate_entropy());
    t_pre_entropy_vec.push_back(t_pre_entropy->calculate_entropy());
  } else if (t_pre_entropy) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      t_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }

    t_pre_entropy->synchronize_entropy_across_nodes();
    t_pre_entropy_vec.push_back(t_pre_entropy->calculate_entropy());
  }

  if (t_pre_entropy && empty_entropy_bins) {
    empty_entropy_frac_vec.push_back(t_pre_entropy->calculate_empty_fraction());
  }
}

void PowerIterator::compute_post_cancellation_entropy(
    std::vector<BankedParticle>& next_gen) {
  if (t_pre_entropy && cancelator) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      p_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      n_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
      t_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }

    p_post_entropy->synchronize_entropy_across_nodes();
    n_post_entropy->synchronize_entropy_across_nodes();
    t_post_entropy->synchronize_entropy_across_nodes();

    p_post_entropy_vec.push_back(p_post_entropy->calculate_entropy());
    n_post_entropy_vec.push_back(n_post_entropy->calculate_entropy());
    t_post_entropy_vec.push_back(t_post_entropy->calculate_entropy());
  }
}

void PowerIterator::compute_pair_dist_sqrd(
    const std::vector<BankedParticle>& next_gen) {
  double Ntot = 0.;
  double r_sqr = 0.;

#ifdef ABEILLE_USE_OMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < next_gen.size(); i++) {
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    Ntot += next_gen[i].wgt;
    for (std::size_t j = 0; j < next_gen.size(); j++) {
      Vector r = next_gen[i].r - next_gen[j].r;

#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
      r_sqr += r.dot(r) * next_gen[i].wgt * next_gen[j].wgt;
    }
  }

  // This is the pair distance squared for the rank.
  r_sqrd /= 2. * Ntot * Ntot;

  // We now get the average pair distance squared across all ranks.
  // This of course leads to a different result with differen parallelization
  // topologies, but is the best that can be done here. This is not a common
  // quantity anyway.
  r_sqr /= static_cast<double>(mpi::size);
  mpi::Allreduce_sum(r_sqr);

  // Save the average pair distance squared for the generation to the variable
  // for this gen in the class, and in the list of all pair distances over the
  // generations.
  r_sqrd = r_sqr;
  r_sqrd_vec.push_back(r_sqrd);
}

void PowerIterator::zero_entropy() {
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

void PowerIterator::premature_kill() {
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
    ngenerations = gen;
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

bool PowerIterator::out_of_time(std::size_t gen) {
  // Get the average time per generation
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

void PowerIterator::check_time(std::size_t gen) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(gen);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    ngenerations = gen;
  }
}

std::shared_ptr<PowerIterator> make_power_iterator(const YAML::Node& sim) {
  // Get the number of particles
  if (!sim["nparticles"] || sim["nparticles"].IsScalar() == false) {
    fatal_error("No nparticles entry in k-eigenvalue simulation.");
  }
  std::size_t nparticles = sim["nparticles"].as<std::size_t>();

  // Get the number of generations
  if (!sim["ngenerations"] || sim["ngenerations"].IsScalar() == false) {
    fatal_error("No ngenerations entry in k-eigenvalue simulation.");
  }
  std::size_t ngenerations = sim["ngenerations"].as<std::size_t>();

  // Get the number of ignored generations
  if (!sim["nignored"] || sim["nignored"].IsScalar() == false) {
    fatal_error("No nignored entry in k-eigenvalue simulation.");
  }
  std::size_t nignored = sim["nignored"].as<std::size_t>();

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
      fatal_error("No sources entry in k-eigenvalue simulation.");
    }
    for (std::size_t s = 0; s < sim["sources"].size(); s++) {
      sources.push_back(make_source(sim["sources"][s]));
    }
  }

  // We now need to make the particle mover, which is the combination of the
  // transport operator and the collision operator.
  std::shared_ptr<IParticleMover> mover =
      make_particle_mover(sim, settings::SimMode::KEFF);

  // Combing
  bool combing = false;
  if (sim["combing"] && sim["combing"].IsScalar()) {
    combing = sim["combing"].as<bool>();
  } else if (sim["combing"]) {
    fatal_error("Invalid combing entry in k-eigenvalue simulation.");
  }

  // Families
  bool families = false;
  if (sim["families"] && sim["families"].IsScalar()) {
    families = sim["families"].as<bool>();
  } else if (sim["families"]) {
    fatal_error("Invalid families entry in k-eigenvalue simulation.");
  }

  // Save source
  bool save_source = false;
  if (sim["save-source"] && sim["save-source"].IsScalar()) {
    families = sim["save-source"].as<bool>();
  } else if (sim["save-source"]) {
    fatal_error("Invalid save-source entry in k-eigenvalue simulation.");
  }

  // Pair Distance Sqrd
  bool pair_dist = false;
  if (sim["pair-distance-sqrd"] && sim["pair-distance-sqrd"].IsScalar()) {
    pair_dist = sim["pair-distance-sqrd"].as<bool>();
  } else if (sim["pair-distance-sqrd"]) {
    fatal_error("Invalid pair-distance-sqrd entry in k-eigenvalue simulation.");
  }

  // Get empty entropy bin fraction
  bool empty_entropy_frac = false;
  if (sim["empty-entropy-bins"] && sim["empty-entropy-bins"].IsScalar()) {
    empty_entropy_frac = sim["empty-entropy-bins"].as<bool>();
  } else if (sim["empty-entropy-bins"]) {
    fatal_error("Invalid empty-entropy-bins entry in k-eigenvalue simulation.");
  }

  // Get entropy bins
  std::optional<Entropy> entropy = std::nullopt;
  if (sim["entropy"] && sim["entropy"].IsMap()) {
    entropy = make_entropy(sim["entropy"]);
  } else if (sim["entropy"]) {
    fatal_error("Invalid entropy entry in k-eigenvalue simulation.");
  }

  // Get cancelator
  std::shared_ptr<Cancelator> cancelator = nullptr;
  if (sim["cancelator"] && sim["cancelator"].IsMap()) {
    cancelator = make_cancelator(sim["cancelator"]);
  } else if (sim["cancelator"]) {
    fatal_error("Invalid cancelator entry in k-eigenvalue simulation.");
  }

  // Make sure our transport operator is compatible with cancellation !
  if (cancelator) cancelator->check_particle_mover_compatibility(mover);

  // Create simulation
  std::shared_ptr<PowerIterator> simptr =
      std::make_shared<PowerIterator>(mover);

  // Set quantities
  simptr->set_nparticles(nparticles);
  simptr->set_ngenerations(ngenerations);
  simptr->set_nignored(nignored);
  for (const auto& src : sources) {
    simptr->add_source(src);
  }
  if (in_source_file.empty() == false) {
    simptr->set_in_source_file(in_source_file);
  }
  simptr->set_combing(combing);
  simptr->set_families(families);
  simptr->set_save_source(save_source);
  simptr->set_pair_distance(pair_dist);
  simptr->set_empty_entropy_bins(empty_entropy_frac);
  if (cancelator) simptr->set_cancelator(cancelator);
  if (entropy) simptr->set_entropy(*entropy);

  return simptr;
}
