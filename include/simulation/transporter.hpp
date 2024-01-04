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
#ifndef TRANSPORTER_H
#define TRANSPORTER_H

#include <geometry/geometry.hpp>
#include <materials/nuclide.hpp>
#include <noise_source/noise_maker.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>
#include <utils/constants.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

#include <optional>

class Transporter {
 public:
  Transporter(std::shared_ptr<Tallies> i_t) : tallies{i_t} {};
  virtual ~Transporter() = default;

  virtual std::vector<BankedParticle> transport(
      std::vector<Particle>& bank, bool noise = false,
      std::vector<BankedParticle>* noise_bank = nullptr,
      const NoiseMaker* noise_maker = nullptr) = 0;

 protected:
  std::shared_ptr<Tallies> tallies;

  struct ThreadLocalScores {
    double k_col_score = 0.;
    double k_abs_score = 0.;
    double k_trk_score = 0.;
    double k_tot_score = 0.;
    double leakage_score = 0.;
    double mig_score = 0.;
  };

  void russian_roulette(Particle& p);

  void collision(Particle& p, MaterialHelper& mat,
                 ThreadLocalScores& thread_scores, bool noise = false,
                 const NoiseMaker* noise_maker = nullptr);

  void branchless_collision(Particle& p, MaterialHelper& mat,
                            ThreadLocalScores& thread_scores);

  void branchless_collision_mat(Particle& p, MaterialHelper& mat,
                                ThreadLocalScores& thread_scores);

  void branchless_collision_iso(Particle& p, MaterialHelper& mat,
                                ThreadLocalScores& thread_scores);

  void branching_collision(Particle& p, MaterialHelper& mat,
                           ThreadLocalScores& thread_scores, bool noise);

  void make_noise_copy(Particle& p, const MicroXSs& microxs) const;

  void do_scatter(Particle& p, const Nuclide& nuclide,
                  const MicroXSs& microxs) const;

  void make_fission_neutrons(Particle& p, const MicroXSs& microxs,
                             const Nuclide& nuclide, bool noise) const;

};  // Transporter

#endif  // MG_TRANSPORTER_H
