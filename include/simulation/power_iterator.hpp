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
#ifndef POWER_ITERATOR_H
#define POWER_ITERATOR_H

#include <simulation/cancelator.hpp>
#include <simulation/simulation.hpp>
#include <unordered_set>
#include <vector>

class PowerIterator : public Simulation {
 public:
  PowerIterator(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
                std::vector<std::shared_ptr<Source>> src)
      : Simulation(i_t, i_tr, src), bank(){};
  PowerIterator(std::shared_ptr<Tallies> i_t, std::shared_ptr<Transporter> i_tr,
                std::vector<std::shared_ptr<Source>> src,
                std::shared_ptr<Cancelator> cncl)
      : Simulation(i_t, i_tr, src), bank(), cancelator(cncl){};
  ~PowerIterator() = default;

  void initialize() override final;

  void run() override final;

  void premature_kill() override final;

 private:
  std::vector<Particle> bank;
  std::shared_ptr<Cancelator> cancelator = nullptr;
  std::unordered_set<uint64_t> families = {};
  int Nnet = 0, Npos = 0, Nneg = 0, Ntot = 0;
  int Wnet = 0, Wpos = 0, Wneg = 0, Wtot = 0;
  int gen = 0;
  double r_sqrd = 0.;

  void check_time(int gen);

  bool out_of_time(int gen);

  void generation_output();

  void write_source(std::vector<Particle>& bank, std::string source_fname);

  void normalize_weights(std::vector<BankedParticle>& next_gen);

  void perform_regional_cancellation(std::vector<BankedParticle>& next_gen);

  // Entropy calculation
  void compute_pre_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void compute_post_cancellation_entropy(std::vector<BankedParticle>& next_gen);
  void zero_entropy();

  // Pair-distance sqrd calculation
  double compute_pair_dist_sqrd(const std::vector<BankedParticle>& next_gen);

  void print_header();

  void sample_source_from_sources();
  void load_source_from_file();

};  // PowerIterator

#endif  // MG_POWER_ITERATOR_H
