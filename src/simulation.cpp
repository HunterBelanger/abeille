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
#include <simulation/simulation.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

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
      p_post_entropy(nullptr),
      n_post_entropy(nullptr),
      t_post_entropy(nullptr) {
  settings::initialize_global_rng();
}

std::vector<Particle> Simulation::sample_sources(int N) {
  // Vector of source weights
  std::vector<double> wgts;
  for (size_t i = 0; i < sources.size(); i++) wgts.push_back(sources[i]->wgt());

  // Generate source particles
  std::vector<Particle> source_particles;
  for (int i = 0; i < N; i++) {
    uint64_t history_id = histories_counter++;
    pcg32 rng(settings::rng_seed);
    uint64_t n_advance = settings::rng_stride * history_id;
    rng.advance(n_advance);
    pcg32 initial_rng = rng;

    int indx = RNG::discrete(rng, wgts);
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
  uint64_t pre_Ntot = bank.size();
  mpi::Allreduce_sum(pre_Ntot);

  // First, we send all particles to the master
  mpi::Gatherv(bank, 0);

  // Now we sort the particles
  if (mpi::rank == 0) {
    std::sort(bank.begin(), bank.end());
  }

  // Now we redistribute particles
  mpi::Scatterv(bank, 0);

  uint64_t post_Ntot = bank.size();
  mpi::Allreduce_sum(post_Ntot);

  if (pre_Ntot != post_Ntot)
    Output::instance()->write("\n post_Ntot != pre_Ntot\n");

  // Make sure each node know how many particles the other has
  nums[mpi::rank] = bank.size();
  for (int n = 0; n < mpi::size; n++) {
    mpi::Bcast<uint64_t>(nums[n], n);
  }
}

void Simulation::particles_to_master(std::vector<BankedParticle>& bank) {
  // First, we send all particles to the master
  mpi::Gatherv(bank, 0);

  // Now we sort the particles
  if (mpi::rank == 0) {
    std::sort(bank.begin(), bank.end());
  } else {
    bank.clear();
  }
}

void Simulation::distribute_particles(std::vector<uint64_t>& nums,
                                      std::vector<BankedParticle>& bank) {
  // Now we redistribute particles
  mpi::Scatterv(bank, 0);

  // Make sure each node know how many particles the other has
  nums[mpi::rank] = bank.size();
  for (int n = 0; n < mpi::size; n++) {
    mpi::Bcast<uint64_t>(nums[n], n);
  }
}
