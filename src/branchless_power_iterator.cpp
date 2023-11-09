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
#include <simulation/branchless_power_iterator.hpp>
#include <utils/error.hpp>
#include <utils/kahan.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

#include <highfive/H5File.hpp>
namespace H5 = HighFive;

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

void BranchlessPowerIterator::initialize() {
  // If not using source file
  if (settings::load_source_file == false) {
    sample_source_from_sources();
  } else {
    // Load source file
    load_source_from_file();
  }

  // Assign History ID as family ID. At beginning, history IDs should start at
  // 1, so this SHOULD yield the desired results.
  for (auto& p : bank) p.set_family_id(p.history_id());
}

void BranchlessPowerIterator::load_source_from_file() {
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
      static_cast<uint64_t>(Nprt / static_cast<uint64_t>(mpi::size));
  uint64_t remainder =
      static_cast<uint64_t>(Nprt % static_cast<uint64_t>(mpi::size));

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

void BranchlessPowerIterator::sample_source_from_sources() {
  Output::instance().write(" Generating source particles...\n");
  // Calculate the base number of particles per node to run
  uint64_t base_particles_per_node =
      static_cast<uint64_t>(settings::nparticles / mpi::size);
  uint64_t remainder = static_cast<uint64_t>(settings::nparticles % mpi::size);

  // Set the base number of particles per node in the node_nparticles vector
  mpi::node_nparticles.resize(static_cast<uint64_t>(mpi::size),
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

void BranchlessPowerIterator::print_header() {
  Output& out = Output::instance();
  out.write("\n");
  std::stringstream output;

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(settings::ngenerations).size());

  if (settings::regional_cancellation && t_pre_entropy) {
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

  if (!settings::regional_cancellation && t_pre_entropy) {
    output << "    Entropy  ";
  }

  if (settings::regional_cancellation && t_pre_entropy) {
    output << "     Positive     Nevative        Total     Positive     "
              "Negative        Total";
  }

  if (settings::regional_cancellation) {
    output << "      Nnet      Ntot      Npos      Nneg      Wnet      Wtot    "
              "  Wpos      Wneg";
  }

  if (settings::pair_distance_sqrd) {
    output << "      <r^2>  ";
  }

  if (settings::families) {
    output << "   N_Families  ";
  }

  if (t_pre_entropy && settings::empty_entropy_bins) {
    output << "   Empty Entropy Frac.  ";
  }

  // Add line underneath
  output << "\n " << std::string(output.str().size() - 3, '-') << "\n";

  out.write(output.str());
}

void BranchlessPowerIterator::generation_output() {
  std::stringstream output;

  double kcol = tallies->kcol();

  if (gen == 1) print_header();

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(settings::ngenerations).size());

  output << " " << std::setw(std::max(n_col_gen, 3)) << std::setfill(' ');
  output << std::right << gen << "   " << std::fixed << std::setprecision(5);
  output << kcol;
  if (gen <= settings::nignored + 1) {
    output << "                      ";
  } else {
    double kcol_avg = tallies->kcol_avg();
    double kcol_err = tallies->kcol_err();
    output << "   " << kcol_avg << " +/- " << kcol_err;
  }
  output << "   ";

  if (!settings::regional_cancellation && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4);
    output << t_pre_entropy_vec.back();
  }

  if (settings::regional_cancellation && t_pre_entropy) {
    output << std::setw(10) << std::scientific << std::setprecision(4)
           << p_pre_entropy_vec.back() << "   ";
    output << n_pre_entropy_vec.back() << "   ";
    output << t_pre_entropy_vec.back() << "   ";

    output << p_post_entropy_vec.back() << "   ";
    output << n_post_entropy_vec.back() << "   ";
    output << t_post_entropy_vec.back() << "   ";
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

  if (settings::pair_distance_sqrd) {
    output << "   " << std::setw(10) << std::scientific << std::setprecision(4)
           << r_sqrd;
  }

  if (settings::families) {
    std::size_t nfamilies = families.size();
    mpi::Reduce_sum(nfamilies, 0);
    families_vec.push_back(nfamilies);
    output << "   " << std::fixed << nfamilies;
  }

  if (t_pre_entropy && settings::empty_entropy_bins) {
    output << "   " << std::setprecision(6) << empty_entropy_frac_vec.back();
  }

  // Add line underneath
  output << "\n";

  Output::instance().write(output.str());
}

void BranchlessPowerIterator::run() {
  Output& out = Output::instance();
  out.write(" Running Branchless k Eigenvalue Problem...\n");
  out.write(" NPARTICLES: " + std::to_string(settings::nparticles) + ", ");
  out.write(" NGENERATIONS: " + std::to_string(settings::ngenerations) + ",");
  out.write(" NIGNORED: " + std::to_string(settings::nignored) + "\n");

  // Zero all entropy bins
  zero_entropy();

  // Start timer
  simulation_timer.reset();
  mpi::synchronize();
  simulation_timer.start();

  // Check for imediate convergence (i.e. we read source from file)
  if (settings::nignored == 0) settings::converged = true;

  for (int g = 1; g <= settings::ngenerations; g++) {
    gen = g;

    if (settings::families) {
      // Before transport, we get the set of all families, so that we can know
      // how many particle families entered the generation.
      families.clear();
      for (auto& p : bank) families.insert(p.family_id());
    }

    std::vector<BankedParticle> next_gen = transporter->transport(bank);
    transported_histories += bank.size();

    if (next_gen.empty()) {
      fatal_error("No fission neutrons were produced.");
    }

    // Do all Pre-Cancelation entropy calculations
    compute_pre_cancellation_entropy(next_gen);

    // Get new keff
    tallies->calc_gen_values();

    mpi::synchronize();

    // Do weight cancelation
    if (settings::regional_cancellation && cancelator) {
      perform_regional_cancellation(next_gen);
    }

    // Calculate net positive and negative weight
    normalize_weights(next_gen);

    // Comb particles
    if (settings::branchless_combing) comb_particles(next_gen);

    // Synchronize particles across nodes. Must be done after combing to
    // maintain reproducibility.
    sync_banks(mpi::node_nparticles, next_gen);

    // Compute pair distance squared
    if (settings::pair_distance_sqrd) compute_pair_dist_sqrd(next_gen);

    // Score the source and gen if passed ignored
    if (settings::converged) {
      tallies->score_source(next_gen, settings::converged);
      tallies->record_generation();
    }

    // Zero tallies for next generation
    tallies->clear_generation();

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
    // recorded, and cancellation has occured. Otherwize source file
    // will be larger than necessary, and wrong number of gens will be
    // in output file based on number averaged for tallies.
    sync_signaled();
    if (signaled) premature_kill();

    // Check if we have enough time to finish the simulation. If not,
    // stop now.
    check_time(g);

    // Zero all entropy bins
    zero_entropy();

    // Once ignored generations are finished, mark as true to start
    // doing tallies
    if (g == settings::nignored) settings::converged = true;
  }

  // Stop timer
  simulation_timer.stop();
  out.write("\n Total Simulation Time: " +
            std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  if (settings::converged) {
    // Write the final results of all estimators
    out.write("\n");
    std::stringstream output;
    output << " Results using " << gen - settings::nignored
           << " active generations:\n";
    output << " -----------------------------------\n";
    output << std::fixed << std::setprecision(6);
    output << " | kcol    = " << tallies->kcol_avg() << " +/- "
           << tallies->kcol_err() << " |\n";
    if (settings::energy_mode == settings::EnergyMode::CE) {
      output << " | kabs    = " << tallies->kabs_avg() << " +/- "
             << tallies->kabs_err() << " |\n";
    }
    if (settings::tracking == settings::TrackingMode::SURFACE_TRACKING) {
      output << " | ktrk    = " << tallies->ktrk_avg() << " +/- "
             << tallies->ktrk_err() << " |\n";
    }
    output << " | leakage = " << tallies->leakage_avg() << " +/- "
           << tallies->leakage_err() << " |\n";
    output << " -----------------------------------\n";

    out.write(output.str());

    // Write saved warnings
    out.write_saved_warnings();

    // Save other outputs
    out.write("\n");
  }

  // Write results to file
  tallies->write_tallies();
  write_entropy_families_etc_to_results();

  // Write source
  write_source(bank);
}

void BranchlessPowerIterator::write_entropy_families_etc_to_results() const {
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  if (families_vec.size() > 0) {
    h5.createDataSet("results/families", families_vec);
  }

  if (r_sqrd_vec.size() > 0) {
    h5.createDataSet("results/pair-dist-sqrd", r_sqrd_vec);
  }

  if (t_pre_entropy_vec.size() > 0 &&
      settings::regional_cancellation == false) {
    h5.createDataSet("results/entropy", t_pre_entropy_vec);
  } else if (t_pre_entropy_vec.size() > 0) {
    h5.createDataSet("results/total-pre-cancel-entropy", t_pre_entropy_vec);
    h5.createDataSet("results/neg-pre-cancel-entropy", n_pre_entropy_vec);
    h5.createDataSet("results/pos-pre-cancel-entropy", p_pre_entropy_vec);

    h5.createDataSet("results/total-post-cancel-entropy", t_post_entropy_vec);
    h5.createDataSet("results/neg-post-cancel-entropy", n_post_entropy_vec);
    h5.createDataSet("results/pos-post-cancel-entropy", p_post_entropy_vec);
  }

  if (settings::empty_entropy_bins) {
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
}

void BranchlessPowerIterator::normalize_weights(
    std::vector<BankedParticle>& next_gen) {
  double W = 0.;
  double W_neg = 0.;
  double W_pos = 0.;
  Nnet = 0;
  Ntot = 0;
  std::tuple<double,double,int,int> tup = kahanBabushkaTuple(next_gen.begin(),next_gen.end());

  W_pos = get<0>(tup);
  W_neg = get<1>(tup);
  Npos = get<2>(tup);
  Nneg = get<3>(tup);
  // Get totals across all nodes
  mpi::Allreduce_sum(W_pos);
  mpi::Allreduce_sum(W_neg);
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
void BranchlessPowerIterator::comb_particles(
    std::vector<BankedParticle>& next_gen) {
  // Get sum of total positive and negative weights using Kahan Summation

  std::tuple<double,double,int,int> tup = kahanBabushkaTuple(next_gen.begin(),next_gen.end());
  double Wpos_node = get<0>(tup);
  double Wneg_node = get<1>(tup);
  

  // Get total weights for all nodes
  std::vector<double> Wpos_each_node(mpi::size, 0.);
  Wpos_each_node[mpi::rank] = Wpos_node;
  mpi::Allreduce_sum(Wpos_each_node);
  const double Wpos = kahanBabushka(Wpos_each_node.begin(), Wpos_each_node.end(), 0.);

  std::vector<double> Wneg_each_node(mpi::size, 0.);
  Wneg_each_node[mpi::rank] = Wneg_node;
  mpi::Allreduce_sum(Wneg_each_node);
  const double Wneg = kahanBabushka(Wneg_each_node.begin(), Wneg_each_node.end(), 0.);

  // Determine how many positive and negative on node and globaly
  const std::size_t Npos_node = static_cast<std::size_t>(std::round(Wpos_node));
  const std::size_t Nneg_node =
      static_cast<std::size_t>(std::round(std::abs(Wneg_node)));

  const std::size_t Npos = static_cast<std::size_t>(std::round(Wpos));
  const std::size_t Nneg = static_cast<std::size_t>(std::round(std::abs(Wneg)));

  // The + 2 is to account for rounding in the ceil operations between global
  // array vs just the node
  std::vector<BankedParticle> new_next_gen;
  new_next_gen.reserve(Npos_node + Nneg_node + 2);

  // Variables for combing particles
  const double avg_pos_wgt = Wpos / static_cast<double>(Npos);
  const double avg_neg_wgt = std::abs(Wneg) / static_cast<double>(Nneg);
  double comb_position_pos = 0.;
  double comb_position_neg = 0.;
  if (mpi::rank == 0) {
    comb_position_pos = RNG::rand(settings::rng) * avg_pos_wgt;
    comb_position_neg = RNG::rand(settings::rng) * avg_neg_wgt;
  }
  mpi::Bcast(comb_position_pos, 0);
  mpi::Bcast(comb_position_neg, 0);

  // Find Comb Position for this Node
  const double U_pos = std::accumulate(Wpos_each_node.begin(),
                                       Wpos_each_node.begin() + mpi::rank, 0.);
  const double beta_pos = std::floor((U_pos - comb_position_pos) / avg_pos_wgt);

  const double U_neg = std::abs(std::accumulate(
      Wneg_each_node.begin(), Wneg_each_node.begin() + mpi::rank, 0.));
  const double beta_neg = std::floor((U_neg - comb_position_neg) / avg_neg_wgt);

  if (mpi::rank != 0) {
    comb_position_pos =
        ((beta_pos + 1.) * avg_pos_wgt) + comb_position_pos - U_pos;
    if (comb_position_pos > avg_pos_wgt) comb_position_pos -= avg_pos_wgt;

    comb_position_neg =
        ((beta_neg + 1.) * avg_neg_wgt) + comb_position_neg - U_neg;
    if (comb_position_neg > avg_neg_wgt) comb_position_neg -= avg_neg_wgt;
  }

  double current_particle_pos = 0.;
  double current_particle_neg = 0.;
  // Comb positive and negative particles at the same time
  for (std::size_t i = 0; i < next_gen.size(); i++) {
    if (next_gen[i].wgt > 0.) {
      current_particle_pos += next_gen[i].wgt;
      while (comb_position_pos < current_particle_pos) {
        new_next_gen.push_back(next_gen[i]);
        new_next_gen.back().wgt = avg_pos_wgt;
        comb_position_pos += avg_pos_wgt;
      }
    } else {
      current_particle_neg -= next_gen[i].wgt;
      while (comb_position_neg < current_particle_neg) {
        new_next_gen.push_back(next_gen[i]);
        new_next_gen.back().wgt = -avg_neg_wgt;
        comb_position_neg += avg_neg_wgt;
      }
    }
  }

  next_gen.swap(new_next_gen);
}

void BranchlessPowerIterator::compute_pre_cancellation_entropy(
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

  if (t_pre_entropy && settings::empty_entropy_bins) {
    empty_entropy_frac_vec.push_back(t_pre_entropy->calculate_empty_fraction());
  }
}

void BranchlessPowerIterator::compute_post_cancellation_entropy(
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

    p_post_entropy_vec.push_back(p_post_entropy->calculate_entropy());
    n_post_entropy_vec.push_back(n_post_entropy->calculate_entropy());
    t_post_entropy_vec.push_back(t_post_entropy->calculate_entropy());
  }
}

void BranchlessPowerIterator::compute_pair_dist_sqrd(
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
  r_sqr /= 2. * Ntot * Ntot;

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

void BranchlessPowerIterator::zero_entropy() {
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

void BranchlessPowerIterator::premature_kill() {
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
    settings::ngenerations = gen;
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

bool BranchlessPowerIterator::out_of_time(int gen) {
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

void BranchlessPowerIterator::check_time(int gen) {
  bool should_stop_now = false;
  if (mpi::rank == 0) should_stop_now = out_of_time(gen);

  mpi::Allreduce_or(should_stop_now);

  if (should_stop_now) {
    // Setting ngenerations to gen will cause us to exit the
    // transport loop as if we finished the simulation normally.
    settings::ngenerations = gen;
  }
}

void BranchlessPowerIterator::perform_regional_cancellation(
    std::vector<BankedParticle>& next_gen) {
  // Only perform cancellation if we are master !!
  std::size_t n_lost_boys = 0;

  // Distribute Particles to Cancellation Bins
  for (auto& p : next_gen) {
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

  if (tmp.size() > 0) next_gen.insert(next_gen.begin(), tmp.begin(), tmp.end());

  // All done ! Clear cancelator for next run
  cancelator->clear();
}
