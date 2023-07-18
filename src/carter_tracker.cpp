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
#include <materials/material.hpp>
#include <materials/material_helper.hpp>
#include <simulation/carter_tracker.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/majorant.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/energy_grid.hpp>

#include <ndarray.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

#ifdef ABEILLE_USE_OMP
#include <omp.h>
#endif

// Currently, carter tracking just finds the majorant xs like is done in
// delta tracking. I am not sure yet how we will implement under-estimations
// of the xs in continuous energy, so for now, we just have this.
CarterTracker::CarterTracker(std::shared_ptr<Tallies> i_t)
    : Transporter(i_t), EGrid(nullptr), Esmp(nullptr) {
  Output::instance().write(" Finding majorant cross sections.\n");
  auto Egrid_Emaj_pair = make_majorant_xs();
  EGrid = std::make_shared<pndl::EnergyGrid>(Egrid_Emaj_pair.first);
  Esmp = std::make_shared<pndl::CrossSection>(Egrid_Emaj_pair.second, EGrid, 0);

  if (settings::energy_mode == settings::EnergyMode::MG) {
    for (uint32_t g = 0; g < settings::ngroups; g++) {
      // Get the energy at the mid-point for the group
      size_t i = g * 2;
      double Eg =
          0.5 * (Egrid_Emaj_pair.first[i] + Egrid_Emaj_pair.first[i + 1]);

      double xs = (*Esmp)(Eg);
      Egrid_Emaj_pair.second[i] = xs * settings::sample_xs_ratio[g];
      Egrid_Emaj_pair.second[i + 1] = xs * settings::sample_xs_ratio[g];
    }

    EGrid = std::make_shared<pndl::EnergyGrid>(Egrid_Emaj_pair.first);
    Esmp =
        std::make_shared<pndl::CrossSection>(Egrid_Emaj_pair.second, EGrid, 0);
  }

  if (mpi::rank == 0) {
    // We now create a temporary array, which will hold the sampling xs info,
    // so we can save it to the output file.
    NDArray<double> smp_xs({2, Egrid_Emaj_pair.first.size()});
    for (std::size_t i = 0; i < Egrid_Emaj_pair.first.size(); i++) {
      smp_xs(0, i) = Egrid_Emaj_pair.first[i];
      smp_xs(1, i) = Egrid_Emaj_pair.second[i];
    }
    auto& h5 = Output::instance().h5();
    auto smp_xs_ds =
        h5.createDataSet<double>("sampling-xs", H5::DataSpace(smp_xs.shape()));
    smp_xs_ds.write_raw(&smp_xs[0]);
  }
}

std::vector<BankedParticle> CarterTracker::transport(
    std::vector<Particle>& bank, bool noise,
    std::vector<BankedParticle>* noise_bank, const NoiseMaker* noise_maker) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel
#endif
  {
    // Thread local storage
    ThreadLocalScores thread_scores;

// Transport all particles in for thread
#ifdef ABEILLE_USE_OMP
#pragma omp for schedule(static)
#endif
    for (size_t n = 0; n < bank.size(); n++) {
      // Particle and its personal tracker
      Particle& p = bank[n];
      Tracker trkr(p.r(), p.u());

      // If we got lost, kill the particle
      if (trkr.is_lost()) {
        std::stringstream mssg;
        mssg << "Particle become lost at " << p.r() << ", ";
        mssg << " u = " << p.u() << ", token = " << trkr.surface_token();
        warning(mssg.str());
        p.kill();
      }
      // Only make helper if we aren't lost, to make sure that material isn't
      // a nullptr
      MaterialHelper mat(trkr.material(), p.E());
      if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);

      // auto bound = trkr.boundary();
      while (p.is_alive()) {
        bool had_collision = false;
        bool crossed_boundary = false;
        auto maj_indx = EGrid->get_lower_index(p.E());
        double Esample = Esmp->evaluate(p.E(), maj_indx) + mat.Ew(p.E(), noise);
        p.set_Esmp(Esample);  // Sampling XS saved for cancellation
        double d_coll = RNG::exponential(p.rng, Esample);
        Boundary bound(INF, -1, BoundaryType::Normal);

        // Try moving the distance to collision, and see if we land in a valid
        // material.
        trkr.move(d_coll);
        trkr.get_current();

        if (trkr.is_lost()) {
          // We got lost. This means we probably flew through a boundary
          // condition. We now go back, and actually find the B.C.
          trkr.set_r(p.r());
          trkr.get_current();
          bound = trkr.get_boundary_condition();
          crossed_boundary = true;
        }

        // Score track length tally for boundary distance.
        // This is here because flux-like tallies are allowed with DT.
        // No other quantity should be scored with a TLE, as an error
        // should have been thrown when building all tallies.
        tallies->score_flight(p, std::min(d_coll, bound.distance), mat,
                              settings::converged);

        if (crossed_boundary) {
          if (bound.boundary_type == BoundaryType::Vacuum) {
            p.kill();
            thread_scores.leakage_score += p.wgt();
            Position r_leak = p.r() + bound.distance * p.u();
            thread_scores.mig_score +=
                p.wgt() * (r_leak - p.r_birth()) * (r_leak - p.r_birth());
          } else if (bound.boundary_type == BoundaryType::Reflective) {
            trkr.do_reflection(p, bound);
            // Check if we are lost
            if (trkr.is_lost()) {
              std::stringstream mssg;
              mssg << "Particle " << p.history_id() << ".";
              mssg << p.secondary_id() << " has become lost.\n";
              mssg << "Previous valid coordinates: r = " << p.previous_r();
              mssg << ", u = " << p.previous_u() << ".\n";
              mssg << "Attempted reflection with surface "
                   << geometry::surfaces[static_cast<std::size_t>(
                                             bound.surface_index)]
                          ->id();
              mssg << " at a distance of " << bound.distance << " cm.\n";
              mssg << "Currently lost at r = " << trkr.r()
                   << ", u = " << trkr.u() << ".";
              fatal_error(mssg.str());
            }
            bound = trkr.get_boundary_condition();
          } else {
            fatal_error("Help me, how did I get here ?");
          }
        } else {
          // Update Position
          p.move(d_coll);
          mat.set_material(trkr.material(), p.E());

          // Get true cross section here
          double Et = mat.Et(p.E(), noise);

          if (Esample >= Et) {
            if (RNG::rand(p.rng) < (Et / Esample)) {
              // Flag real collision
              had_collision = true;
            }
          } else {
            double D = Et / (2. * Et - Esample);
            double F = Et / (D * Esample);

            if ((D - RNG::rand(p.rng)) > 0.) {
              p.set_weight(p.wgt() * F);
              had_collision = true;
            } else {
              p.set_weight(-p.wgt() * F);
            }
          }
        }

        if (p.is_alive() && had_collision) {  // real collision
          collision(p, mat, thread_scores, noise, noise_maker);
          trkr.set_u(p.u());
          p.set_previous_collision_real();
          if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);
        } else if (p.is_alive()) {
          p.set_previous_collision_virtual();
        }

        if (p.is_alive() && std::abs(p.wgt()) >= settings::wgt_split) {
          // Split particle is weight magnitude is too large
          int n_new = static_cast<int>(std::ceil(std::abs(p.wgt())));
          p.split(n_new);
        }

        if (!p.is_alive()) {
          // Attempt a resurection
          p.resurect();

          if (p.is_alive()) {
            trkr.set_r(p.r());
            trkr.set_u(p.u());
            trkr.restart_get_current();
            // Check if we are lost
            if (trkr.is_lost()) {
              std::stringstream mssg;
              mssg << "Particle " << p.history_id() << ".";
              mssg << p.secondary_id() << " has become lost.\n";
              mssg << "Attempted resurection at r = " << trkr.r();
              mssg << ", u = " << trkr.u() << ".";
              fatal_error(mssg.str());
            }
            bound = trkr.get_boundary_condition();
            mat.set_material(trkr.material(), p.E());
            if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);
          } else if (settings::rng_stride_warnings) {
            // History is truly dead.
            // Check if we went past the particle stride.
            uint64_t n_rng_calls = p.number_of_rng_calls();
            if (n_rng_calls > settings::rng_stride) {
              // This isn't really a good thing. We should
              // write a warning.
              std::string mssg = "History " + std::to_string(p.history_id()) +
                                 " overran the RNG stride.";
              warning(mssg);
            }
          }
        }
      }  // While alive
    }    // For all particles

    // Send all thread local scores to tallies instance
    tallies->score_k_col(thread_scores.k_col_score);
    tallies->score_k_abs(thread_scores.k_abs_score);
    tallies->score_k_trk(thread_scores.k_trk_score);
    tallies->score_k_tot(thread_scores.k_tot_score);
    tallies->score_leak(thread_scores.leakage_score);
    tallies->score_mig_area(thread_scores.mig_score);
    thread_scores.k_col_score = 0.;
    thread_scores.k_abs_score = 0.;
    thread_scores.k_trk_score = 0.;
    thread_scores.k_tot_score = 0.;
    thread_scores.leakage_score = 0.;
    thread_scores.mig_score = 0.;
  }  // Parallel

  // Vector to contain all fission daughters for all threads
  std::vector<BankedParticle> fission_neutrons;

  // Empty all particle fission banks into the main one
  for (auto p : bank) {
    p.empty_fission_bank(fission_neutrons);
  }

  if (noise_bank && noise_maker) {
    for (auto& p : bank) {
      p.empty_noise_bank(*noise_bank);
    }
  }

  // Can now clear the old bank
  bank.clear();

  return fission_neutrons;
}
