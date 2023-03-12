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
#include <cmath>
#include <fstream>
#include <geometry/geometry.hpp>
#include <iomanip>
#include <ios>
#include <iostream>
#include <numeric>
#include <simulation/branchless_power_iterator.hpp>
#include <sstream>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>
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
  std::shared_ptr<Output> out = Output::instance();

  // Load source file
  NDArray<double> source = NDArray<double>::load(settings::in_source_file_name);
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
    double x = source[(file_start_loc + i) * 8 + 0];
    double y = source[(file_start_loc + i) * 8 + 1];
    double z = source[(file_start_loc + i) * 8 + 2];
    double ux = source[(file_start_loc + i) * 8 + 3];
    double uy = source[(file_start_loc + i) * 8 + 4];
    double uz = source[(file_start_loc + i) * 8 + 5];
    double E = source[(file_start_loc + i) * 8 + 6];
    double w = source[(file_start_loc + i) * 8 + 7];

    bank.push_back({{x, y, z}, {ux, uy, uz}, E, w, histories_counter++});
    bank.back().initialize_rng(settings::rng_seed, settings::rng_stride);
    tot_wgt += w;
  }
  global_histories_counter = Nprt;

  mpi::Allreduce_sum(tot_wgt);

  out->write(" Total Weight of System: " + std::to_string(std::round(tot_wgt)) +
             "\n");
  settings::nparticles = static_cast<int>(std::round(tot_wgt));
  tallies->set_total_weight(std::round(tot_wgt));
}

void BranchlessPowerIterator::sample_source_from_sources() {
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Generating source particles...\n");
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
  std::shared_ptr<Output> out = Output::instance();
  out->write("\n");
  std::stringstream output;

  // First get the number of columns required to print the max generation number
  int n_col_gen =
      static_cast<int>(std::to_string(settings::ngenerations).size());

  if (settings::regional_cancellation && t_pre_entropy) {
    out->write(
        std::string(static_cast<std::size_t>(std::max(n_col_gen, 3) + 1), ' ') +
        "                                         Pre-Cancelation Entropies    "
        "         Post-Cancelation Entropies\n");
    out->write(
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

  // Add line underneath
  output << "\n " << std::string(output.str().size() - 3, '-') << "\n";

  out->write(output.str());
}

void BranchlessPowerIterator::generation_output() {
  std::shared_ptr<Output> out = Output::instance();
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
    output << t_post_entropy->calculate_entropy();
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

  if (settings::pair_distance_sqrd) {
    output << "   " << std::setw(10) << std::scientific << std::setprecision(4)
           << r_sqrd;
  }

  if (settings::families) {
    output << "   " << std::fixed << families.size();
  }

  // Add line underneath
  output << "\n";

  out->write(output.str());
}

void BranchlessPowerIterator::run() {
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Running Branchless k Eigenvalue Problem...\n");
  out->write(" NPARTICLES: " + std::to_string(settings::nparticles) + ", ");
  out->write(" NGENERATIONS: " + std::to_string(settings::ngenerations) + ",");
  out->write(" NIGNORED: " + std::to_string(settings::nignored) + "\n");

  // Zero all entropy bins
  zero_entropy();

  // Start timer
  simulation_timer.reset();
  mpi::synchronize();
  simulation_timer.start();

  for (int g = 1; g <= settings::ngenerations; g++) {
    gen = g;

    if (settings::families) {
      // Before transport, we get the set of all families, so that we can know
      // how many particle families entered the generation.
      families.clear();
      for (auto& p : bank) families.insert(p.family_id());
    }

    std::vector<BankedParticle> next_gen = transporter->transport(bank);

    if (next_gen.size() == 0) {
      std::string mssg = "No fission neutrons were produced.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Do all Pre-Cancelation entropy calculations
    compute_pre_cancellation_entropy(next_gen);

    // Get new keff
    tallies->calc_gen_values();

    std::sort(next_gen.begin(), next_gen.end());
    mpi::synchronize();

    // Send all particles to master for cancellation and normalization/combing
    particles_to_master(next_gen);

    // Do weight cancelation
    if (settings::regional_cancellation && cancelator && mpi::rank == 0) {
      perform_regional_cancellation(next_gen);
    }

    //// Calculate net positive and negative weight
    //if (mpi::rank == 0) normalize_weights(next_gen);

    //// Comb particles
    //if (settings::branchless_combing && mpi::rank == 0)
    //  comb_particles(next_gen);
    if (mpi::rank == 0) {
      sample_without_replacement(next_gen); 
    }

    // Compute pair distance squared
    if (settings::pair_distance_sqrd && mpi::rank == 0) {
      r_sqrd = compute_pair_dist_sqrd(next_gen);
    }
    mpi::Bcast(r_sqrd, 0);

    // Send particles back to nodes
    mpi::synchronize();
    distribute_particles(mpi::node_nparticles, next_gen);

    // Score the source
    if (settings::converged) {
      tallies->score_source(next_gen, settings::converged);
    }

    // Store gen if passed ignored
    if (settings::converged) tallies->record_generation();

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
  out->write("\n Total Simulation Time: " +
             std::to_string(simulation_timer.elapsed_time()) + " seconds.\n");

  if (settings::converged) {
    // Write the final results of all estimators
    out->write("\n");
    std::stringstream output;
    output << " Results using " << gen - settings::nignored
           << " active generations:\n";
    output << " -----------------------------------\n";
    output << std::fixed << std::setprecision(6);
    output << " | kcol    = " << tallies->kcol_avg() << " +/- "
           << tallies->kcol_err() << " |\n";
    output << " | kabs    = " << tallies->kabs_avg() << " +/- "
           << tallies->kabs_err() << " |\n";
    if (settings::tracking == settings::TrackingMode::SURFACE_TRACKING) {
      output << " | ktrk    = " << tallies->ktrk_avg() << " +/- "
             << tallies->ktrk_err() << " |\n";
    }
    output << " | leakage = " << tallies->leakage_avg() << " +/- "
           << tallies->leakage_err() << " |\n";
    output << " -----------------------------------\n";

    output << "\n Migration Area = " << tallies->mig_area_avg() << " +/- "
           << tallies->mig_area_err() << " \n";
    out->write(output.str());

    // Write saved warnings
    out->write_saved_warnings();

    // Save other outputs
    out->write("\n");

    // Write flux file
    if (settings::converged) {
      tallies->write_tallies();
    }
  }

  // Check to write source
  if (settings::save_source) write_source(bank, settings::source_file_name);
}

void BranchlessPowerIterator::normalize_weights(
    std::vector<BankedParticle>& next_gen) {
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

void BranchlessPowerIterator::comb_particles(
    std::vector<BankedParticle>& next_gen) {
  std::vector<BankedParticle> positive_particles;
  std::vector<BankedParticle> negative_particles;
  positive_particles.reserve(next_gen.size());
  negative_particles.reserve(next_gen.size() / 3);
  double Wpos = 0.;
  double Wneg = 0.;

  // Sort positive and negative particles into different buffers
  for (std::size_t i = 0; i < next_gen.size(); i++) {
    if (next_gen[i].wgt > 0.) {
      Wpos += next_gen[i].wgt;
      positive_particles.push_back(next_gen[i]);
    } else {
      Wneg += next_gen[i].wgt;
      negative_particles.push_back(next_gen[i]);
    }
  }
  next_gen.clear();

  // Determine how many positive and negative
  std::size_t Npos = static_cast<std::size_t>(std::ceil(Wpos));
  std::size_t Nneg = static_cast<std::size_t>(std::ceil(std::abs(Wneg)));
  next_gen.reserve(Npos + Nneg);

  // Comb Positive particles
  std::shuffle(positive_particles.begin(), positive_particles.end(),
               settings::rng);
  double avg_pos_wgt = Wpos / static_cast<double>(Npos);
  double comb_pos = RNG::rand(settings::rng) * avg_pos_wgt;
  double current_particle = 0.;
  for (std::size_t i = 0; i < positive_particles.size(); i++) {
    current_particle += positive_particles[i].wgt;
    while (comb_pos < current_particle) {
      next_gen.push_back(positive_particles[i]);
      next_gen.back().wgt = avg_pos_wgt;
      comb_pos += avg_pos_wgt;
    }
  }

  // Comb Negative Particles
  std::shuffle(negative_particles.begin(), negative_particles.end(),
               settings::rng);
  double avg_neg_wgt = std::abs(Wneg) / static_cast<double>(Npos);
  comb_pos = RNG::rand(settings::rng) * avg_neg_wgt;
  current_particle = 0.;
  for (std::size_t i = 0; i < negative_particles.size(); i++) {
    current_particle -= negative_particles[i].wgt;
    while (comb_pos < current_particle) {
      next_gen.push_back(positive_particles[i]);
      next_gen.back().wgt = -avg_neg_wgt;
      comb_pos += avg_neg_wgt;
    }
  }

  // Reshuffle full buffer
  std::shuffle(next_gen.begin(), next_gen.end(), settings::rng);
}

void BranchlessPowerIterator::sample_without_replacement(std::vector<BankedParticle>& next_gen) {
  std::vector<BankedParticle> new_next_gen;
  new_next_gen.reserve(static_cast<std::size_t>(settings::nparticles));
  
  // The total number of particles we need to sample 
  std::size_t N_to_sample = static_cast<std::size_t>(settings::nparticles);

  if (next_gen.size() < N_to_sample) {
    // If we don't have at least N_to_sample particles in the bank, we first
    // add however many copies of each which are needed. 
    const std::size_t ncopies = static_cast<std::size_t>(std::floor(static_cast<double>(N_to_sample)/static_cast<double>(next_gen.size())));

    for (const auto& p : next_gen) {
      for (std::size_t i = 0; i < ncopies; i++) {
        new_next_gen.push_back(p); 
        new_next_gen.back().wgt = 1.;
        new_next_gen.back().wgt2 = 0.;
      }
    }

    N_to_sample -= next_gen.size() * ncopies;
  }
  
  // First, get all the particle weights
  std::vector<double> weights(next_gen.size(), 0.);
  for (std::size_t i = 0; i < next_gen.size(); i++) {
    weights[i] = next_gen[i].wgt; 
  }
  
  // We now sample all of our particles
  for (std::size_t i = 0; i < N_to_sample; i++) {
    // Sample index
    std::size_t ip = static_cast<std::size_t>(RNG::discrete(settings::rng, weights));

    // Add index to new bank
    new_next_gen.push_back(next_gen[ip]);
    new_next_gen.back().wgt = 1.;

    // Set the weight of the sampled particle to zero, so that it will never
    // be sampled again (i.e. it has not been replaced).
    weights[ip] = 0.;
  }

  next_gen = new_next_gen;
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
  } else if (t_pre_entropy) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      t_pre_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }
    t_pre_entropy->synchronize_entropy_across_nodes();
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
  } else if (t_pre_entropy) {
    for (size_t i = 0; i < next_gen.size(); i++) {
      t_post_entropy->add_point(next_gen[i].r, next_gen[i].wgt);
    }
    t_post_entropy->synchronize_entropy_across_nodes();
  }
}

double BranchlessPowerIterator::compute_pair_dist_sqrd(
    const std::vector<BankedParticle>& next_gen) {
  double Ntot = 0.;
  double r_sqr = 0.;

#pragma omp parallel for
  for (std::size_t i = 0; i < next_gen.size(); i++) {
#pragma omp atomic
    Ntot += next_gen[i].wgt;
    for (std::size_t j = 0; j < next_gen.size(); j++) {
      Vector r = next_gen[i].r - next_gen[j].r;
#pragma omp atomic
      r_sqr += r.dot(r) * next_gen[i].wgt * next_gen[j].wgt;
    }
  }

  r_sqr /= 2. * Ntot * Ntot;

  return r_sqr;
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

void BranchlessPowerIterator::write_source(std::vector<Particle>& bank,
                                           std::string source_fname) {
  // Convert the vector of particles to a vector of BakedParticle
  std::vector<BankedParticle> tmp_bank(bank.size());
  for (std::size_t i = 0; i < bank.size(); i++) {
    tmp_bank[i].r = bank[i].r();
    tmp_bank[i].u = bank[i].u();
    tmp_bank[i].E = bank[i].E();
    tmp_bank[i].wgt = bank[i].wgt();
  }

  // Send all particles to the master process, so that all fission
  // sites can be written to a single npy file.
  mpi::Gatherv(tmp_bank, 0);

  Output::instance()->write("\n Writing source distribution file...\n");

  if (mpi::rank == 0) {
    // Make an NDArray to contain all particles info first
    NDArray<double> source({tmp_bank.size(), 8});

    // Add all particles to the array
    for (std::size_t i = 0; i < tmp_bank.size(); i++) {
      const auto& p = tmp_bank[i];

      source[i * 8 + 0] = p.r.x();
      source[i * 8 + 1] = p.r.y();
      source[i * 8 + 2] = p.r.z();
      source[i * 8 + 3] = p.u.x();
      source[i * 8 + 4] = p.u.y();
      source[i * 8 + 5] = p.u.z();
      source[i * 8 + 6] = p.E;
      source[i * 8 + 7] = p.wgt;
    }

    source.save(source_fname + ".npy");
  }
}

void BranchlessPowerIterator::premature_kill() {
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
    out->write("\n Simulation has been stopped by user. Cleaning up...");
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
  next_gen.insert(next_gen.end(), tmp.begin(), tmp.end());
  mpi::node_nparticles[0] += tmp.size();

  // All done ! Clear cancelator for next run
  cancelator->clear();
}
