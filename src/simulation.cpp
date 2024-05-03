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
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <ndarray.hpp>

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
    // Make an NDArray to contain all particles info first
    NDArray<double> source({tmp_bank.size(), 9});

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
    auto source_dset =
        h5.createDataSet<double>("source", H5::DataSpace(source.shape()));
    source_dset.write_raw(&source[0]);
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
