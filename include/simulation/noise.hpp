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
#ifndef NOISE_H
#define NOISE_H

#include <memory>
#include <simulation/cancelator.hpp>
#include <simulation/noise_maker.hpp>
#include <simulation/simulation.hpp>

class Noise : public Simulation {
 public:
  Noise(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
        std::vector<std::shared_ptr<Source>> srcs, NoiseMaker noise_mkr)
      : Simulation(i_t, i_tr, srcs),
        noise_maker(noise_mkr),
        bank(),
        noise_bank(),
        noise_timer(),
        cancellation_timer(),
        power_iteration_timer(),
        convergence_timer(),
        noise_batch_timer() {}

  Noise(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
        std::vector<std::shared_ptr<Source>> srcs,
        std::shared_ptr<Cancelator> cancel, NoiseMaker noise_mkr)
      : Simulation(i_t, i_tr, srcs),
        noise_maker(noise_mkr),
        bank(),
        noise_bank(),
        noise_timer(),
        cancellation_timer(),
        power_iteration_timer(),
        convergence_timer(),
        noise_batch_timer(),
        cancelator(cancel) {}

  void initialize() override final;
  void run() override final;
  void premature_kill() override final;

 private:
  NoiseMaker noise_maker;
  std::vector<Particle> bank;
  std::vector<BankedParticle> noise_bank;
  Timer noise_timer;
  Timer cancellation_timer;
  Timer power_iteration_timer;
  Timer convergence_timer;
  Timer noise_batch_timer;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  int noise_batch = 0;  // Counter for number of noise batches
  int pi_gen = 0;       // Counter for number of power iteration generations
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;

  // Private helper methods
  void compute_pre_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void compute_post_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void zero_entropy();

  bool out_of_time(int gen);
  void check_time(int gen);

  void normalize_weights(std::vector<BankedParticle>& next_gen);
  void perform_regional_cancellation(std::vector<uint64_t>& nums,
                                     std::vector<BankedParticle>& bank);

  void power_iteration(bool sample_noise);
  void pi_generation_output();

  void noise_simulation();
  void noise_output();

  void print_header() const;

  void sample_source_from_sources();
  void load_source_from_file();
};

#endif
