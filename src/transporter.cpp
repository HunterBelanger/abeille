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
#include <complex>
#include <simulation/transporter.hpp>
#include <sstream>
#include <stdexcept>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

void Transporter::russian_roulette(Particle& p) {
  // Roulette first weight
  if (std::abs(p.wgt()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt()) / settings::wgt_survival);
    if (RNG::rand(p.rng) < P_kill)
      p.set_weight(0.);
    else {
      p.set_weight(std::copysign(settings::wgt_survival, p.wgt()));
    }
  }

  // Roulette second weight
  if (std::abs(p.wgt2()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt2()) / settings::wgt_survival);
    if (RNG::rand(p.rng) < P_kill)
      p.set_weight2(0.);
    else {
      p.set_weight2(std::copysign(settings::wgt_survival, p.wgt2()));
    }
  }

  // Check if we are still alive
  if (p.wgt() == 0. && p.wgt2() == 0.) p.kill();
}

void Transporter::collision(Particle& p, MaterialHelper& mat,
                            ThreadLocalScores& thread_scores, bool noise,
                            const NoiseMaker* noise_maker) {
  // Score flux collision estimator with Sigma_t
  tallies->score_collision(p, mat, settings::converged);
  if (!noise) {
    double k_col_scr = p.wgt() * mat.vEf(p.E()) / mat.Et(p.E(), noise);
    double mig_dist = (p.r() - p.r_birth()).norm();
    double mig_area_scr =
        p.wgt() * mat.Ea(p.E()) / mat.Et(p.E()) * mig_dist * mig_dist;
    thread_scores.k_col_score += k_col_scr;
    thread_scores.mig_score += mig_area_scr;
  }

  // Sample noise source from the noise maker.
  if (noise_maker) {
    noise_maker->sample_noise_source(p, mat, tallies->keff(),
                                     settings::w_noise);
  }

  switch (settings::mode) {
    case settings::SimulationMode::BRANCHLESS_K_EIGENVALUE: {
      if (noise) {
        fatal_error(
            "Cannot perform noise simulations with branchless collisions.");
      }
      branchless_collision(p, mat, thread_scores);
      break;
    }
    default:
      branching_collision(p, mat, thread_scores, noise);
      break;
  }
}

void Transporter::branchless_collision(Particle& p, MaterialHelper& mat,
                                       ThreadLocalScores& thread_scores) {
  if (settings::branchless_material) {
    branchless_collision_mat(p, mat, thread_scores);
  } else {
    branchless_collision_iso(p, mat, thread_scores);
  }
}

void Transporter::branchless_collision_iso(Particle& p, MaterialHelper& mat,
                                           ThreadLocalScores& thread_scores) {
  // Sample a nuclide for the collision
  std::pair<const Nuclide*, MicroXSs> nuclide_info =
      mat.sample_nuclide(p.E(), p.rng, false);
  const Nuclide& nuclide = *nuclide_info.first;
  const MicroXSs& microxs = nuclide_info.second;

  // Score contribution to kabs
  double k_abs_scr =
      p.wgt() * microxs.nu_total * microxs.fission / microxs.total;
  thread_scores.k_abs_score += k_abs_scr;

  // Calcuate scatter probability and weight multiplier
  const double Pscatter = (microxs.total - microxs.absorption) /
                          (microxs.nu_total * microxs.fission +
                           (microxs.total - microxs.absorption));
  const double m = (microxs.nu_total * microxs.fission +
                    (microxs.total - microxs.absorption)) /
                   microxs.total;

  // Roulette
  russian_roulette(p);

  if (p.is_alive()) {
    if (RNG::rand(p.rng) < Pscatter) {
      // Do scatter
      ScatterInfo sinfo = nuclide.sample_scatter(p.E(), p.u(), microxs, p.rng);

      if (sinfo.yield == 0.) {
        // this scatter had a yield of 0... we have no choice but to kill
        // the particle now.
        p.kill();
        return;
      }

      p.set_energy(sinfo.energy);
      p.set_direction(sinfo.direction);
      p.set_weight(p.wgt() * m * sinfo.yield);
      p.set_weight2(p.wgt2() * m * sinfo.yield);

      // Kill if energy is too low
      if (p.E() < settings::min_energy) {
        p.kill();
      }

      // Split particle if weight magnitude is too large
      if (settings::branchless_splitting && p.is_alive() &&
          std::abs(p.wgt()) >= settings::wgt_split) {
        int n_new = static_cast<int>(std::floor(std::abs(p.wgt())));
        p.split(n_new);
      }
    } else {
      // Do fission
      double P_delayed = microxs.nu_delayed / microxs.nu_total;
      FissionInfo finfo = nuclide.sample_fission(
          p.E(), p.u(), microxs.energy_index, P_delayed, p.rng);

      // Construct the new fission particle
      BankedParticle fiss_particle{p.r(),
                                   finfo.direction,
                                   finfo.energy,
                                   p.wgt() * m,
                                   p.wgt2() * m,
                                   p.history_id(),
                                   p.daughter_counter(),
                                   p.family_id(),
                                   p.previous_collision_virtual(),
                                   p.previous_r(),
                                   p.previous_u(),
                                   p.previous_E(),
                                   p.E(),
                                   p.Esmp()};
      p.add_fission_particle(fiss_particle);
      p.kill();
    }
  }
}

void Transporter::branchless_collision_mat(Particle& p, MaterialHelper& mat,
                                           ThreadLocalScores& thread_scores) {
  // Fist calculate the scatter probability and the weight multiplier
  const double Es = mat.Es(p.E());
  const double vEf = mat.vEf(p.E());
  const double Et = mat.Et(p.E());
  const double Pscatter = Es / (vEf + Es);
  const double m = (vEf + Es) / Et;

  // Sample wether we will scatter or fission
  MaterialHelper::BranchlessReaction reaction =
      MaterialHelper::BranchlessReaction::FISSION;
  if (RNG::rand(p.rng) < Pscatter)
    reaction = MaterialHelper::BranchlessReaction::SCATTER;

  // Sample nuclide based on reaction type
  std::pair<const Nuclide*, MicroXSs> nuclide_info =
      mat.sample_branchless_nuclide(p.E(), p.rng, reaction);
  const Nuclide* nuclide = nuclide_info.first;
  MicroXSs microxs = nuclide_info.second;

  // Score kabs
  const double m_i = (microxs.nu_total * microxs.fission +
                      (microxs.total - microxs.absorption)) /
                     microxs.total;
  const double k_abs_scr =
      (m / m_i) * p.wgt() * microxs.nu_total * microxs.fission / microxs.total;
  thread_scores.k_abs_score += k_abs_scr;

  // Roulette
  russian_roulette(p);
  if (p.is_alive() == false) return;

  if (reaction == MaterialHelper::BranchlessReaction::SCATTER) {
    // Do scatter
    ScatterInfo sinfo = nuclide->sample_scatter(p.E(), p.u(), microxs, p.rng);

    if (sinfo.yield == 0.) {
      // this scatter had a yield of 0... we have no choice but to kill
      // the particle now.
      p.kill();
      return;
    }

    p.set_energy(sinfo.energy);
    p.set_direction(sinfo.direction);
    p.set_weight(p.wgt() * m * sinfo.yield);
    p.set_weight2(p.wgt2() * m * sinfo.yield);

    // Kill if energy is too low
    if (p.E() < settings::min_energy) {
      p.kill();
    }

    // Split particle if weight magnitude is too large
    if (settings::branchless_splitting && p.is_alive() &&
        std::abs(p.wgt()) >= settings::wgt_split) {
      int n_new = static_cast<int>(std::ceil(std::abs(p.wgt())));
      p.split(n_new);
    }
  } else {
    // Do fission
    double P_delayed = microxs.nu_delayed / microxs.nu_total;
    FissionInfo finfo = nuclide->sample_fission(
        p.E(), p.u(), microxs.energy_index, P_delayed, p.rng);

    // Construct the new fission particle
    BankedParticle fiss_particle{p.r(),
                                 finfo.direction,
                                 finfo.energy,
                                 p.wgt() * m,
                                 p.wgt2() * m,
                                 p.history_id(),
                                 p.daughter_counter(),
                                 p.family_id(),
                                 p.previous_collision_virtual(),
                                 p.previous_r(),
                                 p.previous_u(),
                                 p.previous_E(),
                                 p.E(),
                                 p.Esmp()};
    p.add_fission_particle(fiss_particle);
    p.kill();
  }
}

void Transporter::branching_collision(Particle& p, MaterialHelper& mat,
                                      ThreadLocalScores& thread_scores,
                                      bool noise) {
  // Sample a nuclide for the collision
  std::pair<const Nuclide*, MicroXSs> nuclide_info =
      mat.sample_nuclide(p.E(), p.rng, noise);
  const Nuclide* nuclide = nuclide_info.first;
  MicroXSs microxs = nuclide_info.second;

  // Score kabs
  if (!noise) {
    double k_abs_scr =
        p.wgt() * microxs.nu_total * microxs.fission / microxs.total;
    thread_scores.k_abs_score += k_abs_scr;
  }

  // Make all fission neutrons
  make_fission_neutrons(p, microxs, *nuclide, noise);

  // Make a copy of the particle when doing noise transport, modifying
  // the weight of the copy as well.
  if (noise) {
    make_noise_copy(p, microxs);
  }

  // Implicit capture
  p.set_weight(p.wgt() * (1. - (microxs.absorption + microxs.noise_copy) /
                                   microxs.total));
  p.set_weight2(p.wgt2() * (1. - (microxs.absorption + microxs.noise_copy) /
                                     microxs.total));

  // Roulette
  russian_roulette(p);

  // Scatter particle
  if (p.is_alive()) {
    // Perform scatter with nuclide
    do_scatter(p, *nuclide, microxs);

    if (p.E() < settings::min_energy) {
      p.kill();
    }
  }
}

void Transporter::do_scatter(Particle& p, const Nuclide& nuclide,
                             const MicroXSs& microxs) const {
  ScatterInfo sinfo = nuclide.sample_scatter(p.E(), p.u(), microxs, p.rng);

  if (sinfo.yield == 0.) {
    // this scatter had a yield of 0... we have no choice but to kill
    // the particle now.
    p.kill();
    return;
  }

  // If yield is an integral quantity, we produce secondary neutrons in the
  // secondary bank, to keep weights from being too high and causing
  // variance issues.
  if (std::floor(sinfo.yield) == sinfo.yield && sinfo.yield != 1.) {
    int n = static_cast<int>(sinfo.yield) - 1;
    for (int i = 0; i < n; i++) {
      // Sample outgoing info
      ScatterInfo ninfo = nuclide.sample_scatter_mt(
          sinfo.mt, p.E(), p.u(), microxs.energy_index, p.rng);

      p.make_secondary(ninfo.direction, ninfo.energy, p.wgt(), p.wgt2());
    }

    // This is set to 1 so that we have the correct weight for the last
    // neutron that we make outside this if block
    sinfo.yield = 1.;
  }

  p.set_direction(sinfo.direction);
  p.set_energy(sinfo.energy);
  p.set_weight(p.wgt() * sinfo.yield);
  p.set_weight2(p.wgt2() * sinfo.yield);
}

void Transporter::make_noise_copy(Particle& p, const MicroXSs& microxs) const {
  std::complex<double> weight_copy{p.wgt(), p.wgt2()};
  if (microxs.noise_copy / microxs.total + RNG::rand(p.rng) >= 1.) {
    std::complex<double> yield{1., -1. / settings::eta};
    weight_copy *= yield;
    p.make_secondary(p.u(), p.E(), weight_copy.real(), weight_copy.imag());
  }
}

void Transporter::make_fission_neutrons(Particle& p, const MicroXSs& microxs,
                                        const Nuclide& nuclide,
                                        bool noise) const {
  double k_abs_scr =
      p.wgt() * microxs.nu_total * microxs.fission / microxs.total;

  // Determine the number of neutrons to produce
  int n_new = 0;
  switch (settings::mode) {
    case settings::SimulationMode::K_EIGENVALUE:
      // In k-eivenvalue simulations, we normalize particle production
      // by keff so that it stays more or less constant.
      n_new = static_cast<int>(
          std::floor(std::abs(k_abs_scr) / tallies->kcol() + RNG::rand(p.rng)));
      break;

    case settings::SimulationMode::FIXED_SOURCE:
      // In fixed-source problems, we don't normalize the number of
      // particles. Problem MUST BE SUB CRITICAL !
      n_new =
          static_cast<int>(std::floor(std::abs(k_abs_scr) + RNG::rand(p.rng)));
      break;

    case settings::SimulationMode::MODIFIED_FIXED_SOURCE:
      // In fixed-source problems, we don't normalize the number of
      // particles. Problem MUST BE SUB CRITICAL !
      n_new =
          static_cast<int>(std::floor(std::abs(k_abs_scr) + RNG::rand(p.rng)));
      break;

    case settings::SimulationMode::NOISE:
      if (!noise) {
        // This is a power iteration generation. Run as with K_EIGENVALUE
        n_new = static_cast<int>(std::floor(
            std::abs(k_abs_scr) / tallies->kcol() + RNG::rand(p.rng)));
      } else {
        // When transporting noise particles, we don't normalize the number of
        // generated particles by the weight ! We also scale by kcol from the
        // previous power-iteration generation, to make the problem "critical".
        n_new =
            static_cast<int>(std::floor((microxs.nu_total * microxs.fission /
                                         (microxs.total * tallies->keff())) +
                                        RNG::rand(p.rng)));
      }
      break;

    case settings::SimulationMode::BRANCHLESS_K_EIGENVALUE:
      fatal_error("Should never get here.");
      break;
  }

  // Probability of a fission neutron being a delayed neutron.
  double P_delayed = microxs.nu_delayed / microxs.nu_total;

  for (int i = 0; i < n_new; i++) {
    auto finfo = nuclide.sample_fission(p.E(), p.u(), microxs.energy_index,
                                        P_delayed, p.rng);

    // If we are doing transport of non-noise particles, the weight of fission
    // particles should be +/- 1, depending on sign of parent.
    double wgt = p.wgt() > 0. ? 1. : -1.;
    double wgt2 = 0.;

    if (noise) {
      // If parent was a noise particle, we start with their weight
      wgt = p.wgt();
      wgt2 = p.wgt2();

      if (finfo.delayed) {
        // If the sampled fission neutron was delayed, we need to apply a
        // complex yield based on the precursor decay constant.
        std::complex<double> wgt_cmpx{wgt, wgt2};
        double lambda = finfo.precursor_decay_constant;
        double denom =
            (lambda * lambda) + (settings::w_noise * settings::w_noise);
        std::complex<double> mult{lambda * lambda / denom,
                                  -lambda * settings::w_noise / denom};
        wgt_cmpx *= mult;
        wgt = wgt_cmpx.real();
        wgt2 = wgt_cmpx.imag();
      }
    }

    // Construct the new fission particle
    BankedParticle fiss_particle{p.r(),
                                 finfo.direction,
                                 finfo.energy,
                                 wgt,
                                 wgt2,
                                 p.history_id(),
                                 p.daughter_counter(),
                                 p.family_id(),
                                 p.previous_collision_virtual(),
                                 p.previous_r(),
                                 p.previous_u(),
                                 p.previous_E(),
                                 p.E(),
                                 p.Esmp()};

    // Save the new fission particle in appropriate manner,
    // depending on the simulation mode.
    switch (settings::mode) {
      case settings::SimulationMode::FIXED_SOURCE:
        p.make_secondary(fiss_particle.u, fiss_particle.E, fiss_particle.wgt,
                         fiss_particle.wgt2);
        break;

      case settings::SimulationMode::MODIFIED_FIXED_SOURCE:
        p.add_fission_particle(fiss_particle);
        break;

      case settings::SimulationMode::K_EIGENVALUE:
        p.add_fission_particle(fiss_particle);
        break;

      case settings::SimulationMode::NOISE:
        if (!noise || settings::inner_generations) {
          p.add_fission_particle(fiss_particle);
        } else {
          p.make_secondary(fiss_particle.u, fiss_particle.E, fiss_particle.wgt,
                           fiss_particle.wgt2);
        }
        break;

      case settings::SimulationMode::BRANCHLESS_K_EIGENVALUE:
        fatal_error("Should never get here.");
        break;
    }
  }
}
