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
#include <simulation/fixed_source.hpp>
#include <simulation/noise.hpp>
#include <simulation/power_iterator.hpp>
#include <simulation/simulation.hpp>
#include <utils/error.hpp>
#include <utils/kahan.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <xtensor/xtensor.hpp>

#include <functional>

Simulation::Simulation(std::shared_ptr<IParticleMover> i_pm)
    : particle_mover{i_pm}, simulation_timer() {
  settings::initialize_global_rng();
}

void Simulation::sync_signaled() {
  mpi::Allreduce_or(signaled);
  mpi::Allreduce_or(terminate);
}

std::vector<Particle> Simulation::sample_sources(
    const std::vector<std::shared_ptr<Source>>& sources, std::size_t N) {
  // Vector of source weights
  std::vector<double> wgts;
  for (size_t i = 0; i < sources.size(); i++) wgts.push_back(sources[i]->wgt());

  // Generate source particles
  std::vector<Particle> source_particles;
  for (std::size_t i = 0; i < N; i++) {
    uint64_t history_id = histories_counter++;
    RNG rng(settings::rng_seed);
    uint64_t n_advance = settings::rng_stride * history_id;
    rng.advance(n_advance);
    RNG initial_rng = rng;

    std::size_t indx = rng.discrete(wgts);
    source_particles.push_back(sources[indx]->generate_particle(rng));
    source_particles.back().set_history_id(history_id);
    source_particles.back().rng = rng;
    source_particles.back().set_initial_rng(initial_rng);
  }

  return source_particles;
}

void Simulation::sync_banks(std::vector<uint64_t>& nums,
                            std::vector<BankedParticle>& bank) {
  // Make sure all ranks know how many particles all other nodes have
  std::fill(nums.begin(), nums.end(), 0);
  nums[static_cast<std::size_t>(mpi::rank)] = bank.size();
  mpi::Allreduce_sum(nums);

  // These are constants used for determining how many particles should be on
  // each rank. Ntot is the total number of particles across all ranks, base is
  // the minimum number of particles that all ranks should have. Remainder is
  // the number of particles which are left if each rank has base particles. If
  // rank < remainder, then that rank gets one extra particle.
  const uint64_t Ntot =
      std::accumulate(nums.begin(), nums.end(), static_cast<uint64_t>(0));
  const uint64_t base = Ntot / static_cast<uint64_t>(mpi::size);
  const uint64_t remainder = Ntot - (static_cast<uint64_t>(mpi::size) * base);

  // This lambda sends nts particles to the LEFT from rank R to rank R-1.
  // It also updates nums for each rank.
  auto SendLeft = [&bank, &nums](const int R, const int nts) {
    const uint64_t unts = static_cast<uint64_t>(nts);
    // nts is Number to send
    std::vector<BankedParticle> sendBank;
    if (mpi::rank == R) {
      sendBank = {bank.begin(), bank.begin() + nts};
      bank.erase(bank.begin(), bank.begin() + nts);
      mpi::Send(std::span<BankedParticle>(sendBank.begin(), sendBank.end()),
                R - 1);
    } else if (mpi::rank == (R - 1)) {
      sendBank.resize(unts);
      mpi::Recv(std::span<BankedParticle>(sendBank.begin(), sendBank.end()), R);
      bank.insert(bank.end(), sendBank.begin(), sendBank.end());
    }
    nums[static_cast<std::size_t>(R)] -= unts;
    nums[static_cast<std::size_t>(R - 1)] += unts;
  };

  // This lambda sends nts particles to the RIGHT from rank L to rank L+1.
  // It also updates nums for each rank.
  auto SendRight = [&bank, &nums](const int L, const int nts) {
    const uint64_t unts = static_cast<uint64_t>(nts);
    // nts is Number to send
    std::vector<BankedParticle> sendBank;
    if (mpi::rank == L) {
      sendBank = {bank.end() - nts, bank.end()};
      bank.erase(bank.end() - nts, bank.end());
      mpi::Send(std::span<BankedParticle>(sendBank.begin(), sendBank.end()),
                L + 1);
    } else if (mpi::rank == (L + 1)) {
      sendBank.resize(unts);
      mpi::Recv(std::span<BankedParticle>(sendBank.begin(), sendBank.end()), L);
      bank.insert(bank.begin(), sendBank.begin(), sendBank.end());
    }
    nums[static_cast<std::size_t>(L)] -= unts;
    nums[static_cast<std::size_t>(L + 1)] += unts;
  };

  std::function<void(int, int)> SendParticles =
      [&nums, base, remainder, SendRight, SendLeft, &SendParticles](
          const int i, const int ntr_left) {
        // Number of particles on ranks i and i+1
        int nums_i = static_cast<int>(nums[static_cast<std::size_t>(i)]);
        int nums_i1 = static_cast<int>(nums[static_cast<std::size_t>(i + 1)]);

        // ntr is of particles that rank i needs to receive. This includes the
        // particles i needs, in addition to any particles previous ranks might
        // need (ntr_left).
        int ntr = ntr_left + static_cast<int>(base) +
                  (i < static_cast<int>(remainder) ? 1 : 0) - nums_i;

        if (ntr > 0 && ntr > nums_i1) {
          SendParticles(i + 1, ntr);

          // Update Nums after recursive call
          nums_i = static_cast<int>(nums[static_cast<std::size_t>(i)]);
          nums_i1 = static_cast<int>(nums[static_cast<std::size_t>(i + 1)]);
          ntr = ntr_left + static_cast<int>(base) +
                (i < static_cast<int>(remainder) ? 1 : 0) - nums_i;
        }

        if (ntr > 0) {
          // Needs particles ! Get ntr particles from the next rank.
          SendLeft(i + 1, ntr);
        } else if (ntr < 0) {
          // Too many particles ! Send abs(ntr) particles to next rank.
          SendRight(i, -ntr);
        }
      };

  // Call SendParticles for all but the last rank. Since we move from rank 0 up
  // to the last rank, if all ranks but the last have the correct number of
  // particles, then the last rank must also have the correct number.
  for (int i = 0; i < mpi::size - 1; i++) {
    SendParticles(i, 0);
  }
}

void Simulation::normalize_weights(std::size_t nparticles,
                                   std::vector<BankedParticle>& next_gen,
                                   int& Npos, int& Nneg, int& Nnet, int& Ntot,
                                   int& Wpos, int& Wneg, int& Wnet, int& Wtot) {
  double W = 0.;
  double W_neg = 0.;
  double W_pos = 0.;

  std::tie(W_pos, W_neg, Npos, Nneg) =
      kahan_bank(next_gen.begin(), next_gen.end());
  W_neg = std::abs(W_neg);

  // Get totals across all nodes
  mpi::Allreduce_sum(W_pos);
  mpi::Allreduce_sum(W_neg);
  mpi::Allreduce_sum(Npos);
  mpi::Allreduce_sum(Nneg);

  W = W_pos - W_neg;
  Ntot = Npos + Nneg;
  Nnet = Npos - Nneg;

  // Re-Normalize particle weights
  double w_per_part = static_cast<double>(nparticles) / W;
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

void Simulation::comb_particles(std::vector<BankedParticle>& next_gen,
                                int& Npos, int& Nneg, int& Nnet, int& Ntot) {
  double Wpos_node = 0.;
  double Wneg_node = 0.;
  // Get sum of total positive and negative weights using Kahan Summation
  std::tie(Wpos_node, Wneg_node, std::ignore, std::ignore) =
      kahan_bank(next_gen.begin(), next_gen.end());

  // Get total weights for all nodes
  std::vector<double> Wpos_each_node(static_cast<std::size_t>(mpi::size), 0.);
  Wpos_each_node[static_cast<std::size_t>(mpi::rank)] = Wpos_node;
  mpi::Allreduce_sum(Wpos_each_node);
  const double Wpos = kahan(Wpos_each_node.begin(), Wpos_each_node.end(), 0.);

  std::vector<double> Wneg_each_node(static_cast<std::size_t>(mpi::size), 0.);
  Wneg_each_node[static_cast<std::size_t>(mpi::rank)] = Wneg_node;
  mpi::Allreduce_sum(Wneg_each_node);
  const double Wneg = kahan(Wneg_each_node.begin(), Wneg_each_node.end(), 0.);

  // Determine how many positive and negative particles we want to have after
  // combing on this node and globaly as well
  const std::size_t Npos_node = static_cast<std::size_t>(std::round(Wpos_node));
  const std::size_t Nneg_node =
      static_cast<std::size_t>(std::round(std::abs(Wneg_node)));

  const std::size_t Npos_st = static_cast<std::size_t>(std::round(Wpos));
  const std::size_t Nneg_st =
      static_cast<std::size_t>(std::round(std::abs(Wneg)));

  // Update particle numbers due to combing
  Npos = static_cast<int>(Npos_st);
  Nneg = static_cast<int>(Nneg_st);
  Ntot = Npos + Nneg;
  Nnet = Npos - Nneg;

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
    // The initial random offset of the comb is sampled on the master node.
    comb_position_pos = settings::rng() * avg_pos_wgt;
    comb_position_neg = settings::rng() * avg_neg_wgt;
  }
  mpi::Bcast(comb_position_pos, 0);
  mpi::Bcast(comb_position_neg, 0);

  // Find Comb Position for this Node
  const double U_pos =
      kahan(Wpos_each_node.begin(), Wpos_each_node.begin() + mpi::rank, 0.);
  const double beta_pos = std::floor((U_pos - comb_position_pos) / avg_pos_wgt);

  const double U_neg = std::abs(
      kahan(Wneg_each_node.begin(), Wneg_each_node.begin() + mpi::rank, 0.));
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
  // Comb positive and negative particles at the same time, to ensure that
  // particle orders aren't changed for different parallelization schemes.
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

void Simulation::perform_regional_cancellation(
    std::shared_ptr<Cancelator>& cancelator,
    std::vector<BankedParticle>& next_gen) {
  if (cancelator == nullptr) return;

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
  cancelator->perform_cancellation();

  // All particles which were placed into a cancellation bin from next_gen
  // now have modified weights.
  // Now we can get the uniform particles
  auto tmp = cancelator->get_new_particles(settings::rng);

  if (tmp.size() > 0) next_gen.insert(next_gen.begin(), tmp.begin(), tmp.end());

  // All done ! Clear cancelator for next run
  cancelator->clear();
}

void Simulation::write_source(std::vector<Particle>& bank) const {
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

  if (mpi::rank == 0) {
    // Make an array to contain all particle's info first
    xt::xtensor<double, 2> source;
    source.resize({tmp_bank.size(), 9});
    source.fill(0.);

    // Add all particles to the array
    for (std::size_t i = 0; i < tmp_bank.size(); i++) {
      const auto& p = tmp_bank[i];

      source[i * 9 + 0] = p.r.x();
      source[i * 9 + 1] = p.r.y();
      source[i * 9 + 2] = p.r.z();
      source[i * 9 + 3] = p.u.x();
      source[i * 9 + 4] = p.u.y();
      source[i * 9 + 5] = p.u.z();
      source[i * 9 + 6] = p.E;
      source[i * 9 + 7] = p.wgt;
      source[i * 9 + 8] = p.wgt2;
    }

    auto& h5 = Output::instance().h5();
    const std::vector<std::size_t> shape{tmp_bank.size(), 9};
    auto source_dset = h5.createDataSet<double>("source", H5::DataSpace(shape));
    source_dset.write_raw(source.data());
  }
}

std::shared_ptr<Simulation> make_simulation(const YAML::Node& input) {
  // Get the simulation mode
  if (!input["simulation"]) {
    fatal_error("No simulation entry provided in input.");
  } else if (input["simulation"] && input["simulation"].IsMap() == false) {
    fatal_error("Invalid simulation entry provided in input.");
  }
  const YAML::Node& sim = input["simulation"];

  // First, get the string identifying the simulation mode
  if (!sim["mode"]) {
    fatal_error("No mode entry in simulation definition.");
  } else if (sim["mode"].IsScalar() == false) {
    fatal_error("Invalid mode entry in simulation.");
  }
  std::string mode = sim["mode"].as<std::string>();

  std::shared_ptr<Simulation> simulation = nullptr;
  if (mode == "k-eigenvalue") {
    simulation = make_power_iterator(sim);
  } else if (mode == "fixed-source") {
    simulation = make_fixed_source(sim);
  } else if (mode == "modified-fixed-source") {
    fatal_error("Simulation mode modified-fixed-source not yet implemented.");
  } else if (mode == "noise") {
    simulation = make_noise_simulator(sim);
  } else {
    fatal_error("Unknown simulation mode " + mode + ".");
  }

  return simulation;
}
