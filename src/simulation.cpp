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
#include <simulation/simulation.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <functional>

#include <ndarray.hpp>

Simulation::Simulation(std::shared_ptr<Tallies> i_t,
                       std::shared_ptr<Transporter> i_tr,
                       std::vector<std::shared_ptr<Source>> srcs)
    : tallies{i_t},
      transporter{i_tr},
      sources{srcs},
      simulation_timer(),
      p_pre_entropy(nullptr),
      n_pre_entropy(nullptr),
      t_pre_entropy(nullptr),
      p_pre_entropy_vec(),
      n_pre_entropy_vec(),
      t_pre_entropy_vec(),
      p_post_entropy(nullptr),
      n_post_entropy(nullptr),
      t_post_entropy(nullptr),
      p_post_entropy_vec(),
      n_post_entropy_vec(),
      t_post_entropy_vec(),
      empty_entropy_frac_vec() {
  settings::initialize_global_rng();
}

std::vector<Particle> Simulation::sample_sources(std::size_t N) {
  // Vector of source weights
  std::vector<double> wgts;
  for (size_t i = 0; i < sources.size(); i++) wgts.push_back(sources[i]->wgt());

  // Generate source particles
  std::vector<Particle> source_particles;
  for (std::size_t i = 0; i < N; i++) {
    uint64_t history_id = histories_counter++;
    pcg32 rng(settings::rng_seed);
    uint64_t n_advance = settings::rng_stride * history_id;
    rng.advance(n_advance);
    pcg32 initial_rng = rng;

    std::size_t indx = static_cast<std::size_t>(RNG::discrete(rng, wgts));
    source_particles.push_back(sources[indx]->generate_particle(rng));
    source_particles.back().set_history_id(history_id);
    source_particles.back().rng = rng;
    source_particles.back().set_initial_rng(initial_rng);
  }

  return source_particles;
}

void Simulation::sync_signaled() {
  mpi::Allreduce_or(signaled);
  mpi::Allreduce_or(terminate);
}

void Simulation::sync_banks(std::vector<uint64_t>& nums,
                            std::vector<BankedParticle>& bank) {

    nums[static_cast<std::size_t>(mpi::rank)] = bank.size();
    for (int n = 0; n < mpi::size; n++) {
    mpi::Bcast<uint64_t>(nums[static_cast<std::size_t>(n)], n);
  }

  uint64_t Ntot = std::accumulate(nums.begin(),nums.end(),0);
  //std::cout << " NTOT IS : " << Ntot;
  //std::cout << mpi::rank << " has " << bank.size() << " particles.\n";
 const int base = static_cast<int>(Ntot) / mpi::size;
 //std::cout << " BASE IS " << base << "\n";
  const int remainder = static_cast<int>(Ntot) - (mpi::size * base);
 //std::cout << " remainder IS " << remainder << "\n";
  auto SendLeft = [&bank, & nums,base](int R, int nts) {
    std::vector<BankedParticle> sendBank;
 	if (mpi::rank == R) {
    sendBank = {bank.begin(), bank.begin() + nts};
    bank.erase(bank.begin(), bank.begin() + nts);
    mpi::Send(sendBank,R-1);
	} else if (mpi::rank == (R-1)) {
    mpi::Recv(sendBank,R);
    bank.insert(bank.end(),sendBank.begin(),sendBank.end());
	}
  };
  auto SendRight = [&bank, &nums,base](int L, int nts) {
  std::vector<BankedParticle> sendBank;
 	if (mpi::rank == L) {
	  sendBank = {bank.begin() + base, bank.begin() + base + nts};
    bank.erase(bank.begin() + base, bank.begin() + base + nts);
    mpi::Send(sendBank,L+1);
	} else if (mpi::rank == (L+1)) {
	  mpi::Recv(sendBank,L);
    bank.insert(bank.begin(),sendBank.begin(),sendBank.end());
	}
  };
  std::function<void(int,int)> DoMagic = [&nums, SendRight, SendLeft, base, remainder, &DoMagic](int i, int ntr_left){

	const int nums_i = static_cast<int>(nums[i]);
	const int nums_i1 = static_cast<int>(nums[i+1]);
	const int ntr = ntr_left + base + (i < remainder ? 1 : 0) - nums_i; 
  //  std::cout << "ntr calculation: " << ntr_left <<  " + " << base << " + " << remainder << " - " << nums_i << "\n";

  //std::cout << "base: " << base;
  //for(int i = 0 ; i  < mpi::size ; i ++)
   // if(mpi::rank == i)
     // std::cout << i << " size: " << nums_i << " ntr is "<< ntr <<  "\n";
 	if (ntr > 0 && ntr > nums_i1) {
    //std::cout << " am i here\n";
	  DoMagic(i+1, ntr);
	} else if (ntr > 0 && ntr <= nums_i1)  {
		SendLeft(i+1, ntr);	
    nums[i+1] -= ntr;
    nums[i] += ntr;
	} 
  //When in the recursive case, if for example the right node does not have enough elements and you went to the one after
  // when the nodes get sent around your ntr < 0 but you want to send left if the right one has enough elements
 /* else if (ntr < 0 && nums_i1 >= base)
  {
    SendLeft(i, -ntr);	
    nums[i] -= -ntr;
    nums[i-1] +=  -ntr;
  }*/
  else if (ntr < 0) {
		SendRight(i, -ntr);	
    nums[i] -= -ntr;
    nums[i+1] += -ntr;
	}
  };

  for (int i = 0; i < mpi::size - 1; i++)
	  DoMagic(i,0);

   // std::cout << mpi::rank << " has " << bank.size() << " particles.\n";
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
