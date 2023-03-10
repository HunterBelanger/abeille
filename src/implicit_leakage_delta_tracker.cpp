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
#ifdef _OPENMP
#include <omp.h>
#endif

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_tsl_reaction.hpp>
#include <materials/material.hpp>
#include <materials/material_helper.hpp>
#include <materials/ce_nuclide.hpp>
#include <simulation/implicit_leakage_delta_tracker.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <boost/unordered/unordered_flat_map.hpp>

#include <algorithm>
#include <iomanip>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

ImplicitLeakageDeltaTracker::ImplicitLeakageDeltaTracker(
    std::shared_ptr<Tallies> i_t)
    : Transporter(i_t), EGrid(nullptr), Emaj(nullptr) {
  Output::instance()->write(" Finding majorant cross sections.\n");
  // Must first create a unionized energy grid. How this is done depends on
  // whether or not we are in continuous energy or multi-group mode.
  if (settings::energy_mode == settings::EnergyMode::CE) {
    // First, we must construct a unionized energy grid, for all nuclides in the
    // problem. We do this by iterating through all nuclides. In CE mode, we
    // should be able to safely case a Nuclide pointer to a CENuclide pointer.
    std::vector<double> union_energy_grid;

    // Get a map of all URR random values set to 1, for finding the majorant
    boost::unordered_flat_map<uint32_t, std::optional<double>> urr_rands;
    for (const auto& za : zaids_with_urr) urr_rands[za] = 1.;

    // Iterate through all nuclides
    for (const auto& id_nucld_pair : nuclides) {
      Nuclide* nuc = id_nucld_pair.second.get();
      CENuclide* cenuc = static_cast<CENuclide*>(nuc);

      // Get the nuclide's energy grid
      std::vector<double> cenuc_grid = cenuc->cedata()->energy_grid().grid();
      if (cenuc->tsl()) {
        // If we have a TSL, we need to change cenuc_grid to add those points.
        // First, get ref to IncoherentInelastic
        const pndl::STTSLReaction& origII =
            cenuc->tsl()->incoherent_inelastic();
        const pndl::STIncoherentInelastic& II =
            reinterpret_cast<const pndl::STIncoherentInelastic&>(origII);
        const pndl::STTSLReaction& origCE = cenuc->tsl()->coherent_elastic();
        const pndl::STCoherentElastic& CE =
            reinterpret_cast<const pndl::STCoherentElastic&>(origCE);
        std::size_t NEtsl = II.xs().x().size() + (2 * CE.bragg_edges().size());
        cenuc_grid.reserve(cenuc_grid.size() + NEtsl);

        // First, add all II points
        for (std::size_t i = 0; i < II.xs().x().size(); i++) {
          cenuc_grid.push_back(II.xs().x()[i]);
        }

        // Now add all CE bragg edges
        for (std::size_t i = 0; i < CE.bragg_edges().size(); i++) {
          cenuc_grid.push_back(std::nextafter(CE.bragg_edges()[i], 0.));
          cenuc_grid.push_back(CE.bragg_edges()[i]);
        }

        // Now sort the grid
        std::sort(cenuc_grid.begin(), cenuc_grid.end());
      }

      // Now we perform the union operation
      std::vector<double> tmp;
      std::set_union(union_energy_grid.begin(), union_energy_grid.end(),
                     cenuc_grid.begin(), cenuc_grid.end(),
                     std::back_inserter(tmp));
      union_energy_grid = tmp;

      // We now throw the grid into a set temporarily, to remove the duplicates
      std::set<double> tmp_union_energy_grid(union_energy_grid.begin(),
                                             union_energy_grid.end());
      union_energy_grid.assign(tmp_union_energy_grid.begin(),
                               tmp_union_energy_grid.end());
    }

    // With the energy grid unionized, we can go about finding majorants !
    std::vector<double> maj_xs(union_energy_grid.size(), 0.);
    for (const auto& material : materials) {
      // Make a MaterialHelper, to evaluate the material XS for us
      MaterialHelper mat(material.second.get(), union_energy_grid.front());
      mat.set_urr_rand_vals(urr_rands);

      // Now iterate through all energy points. If xs is larger, update value
      for (std::size_t i = 0; i < union_energy_grid.size(); i++) {
        const double Ei = union_energy_grid[i];
        const double xsi = mat.Et(Ei);

        if (xsi > maj_xs[i]) maj_xs[i] = xsi;
      }
    }

    // With the majorant obtained, we can now construct the majorant xs
    EGrid = std::make_shared<pndl::EnergyGrid>(union_energy_grid);
    Emaj = std::make_shared<pndl::CrossSection>(maj_xs, EGrid, 0);
  } else {
    // We are in multi-group mode. Here, the energy-bounds are kept in the
    // settings, so we can construct something with that

    std::vector<double> egrid;
    egrid.push_back(settings::energy_bounds[0]);
    if (egrid.front() == 0.) egrid.front() = 1.E-11;

    for (size_t i = 1; i < settings::energy_bounds.size() - 1; i++) {
      egrid.push_back(settings::energy_bounds[i]);
      egrid.push_back(settings::energy_bounds[i]);
    }

    egrid.push_back(settings::energy_bounds.back());

    // This now has created the vector egrid which will look something like this
    // for the case of 5 energy groups.
    // [0., 1.,   1., 2.,   2., 3.,   3., 4.,   4., 5.]
    // This works, because the energy of multi-group particles should always be
    // inbetween the bounds for the group.

    // Now we need to make a vector which will contian the majorant cross
    // cross sections for each group.
    std::vector<double> maj_xs(egrid.size(), 0.);

    // We loop through materials
    for (const auto& material : materials) {
      // Then we loop through energies
      MaterialHelper mat(material.second.get(), 1.);

      for (uint32_t g = 0; g < settings::ngroups; g++) {
        // Get the energy at the mid-point for the group
        size_t i = g * 2;
        double Eg = 0.5 * (egrid[i] + egrid[i + 1]);

        double xs = mat.Et(Eg);
        if (xs > maj_xs[i]) {
          maj_xs[i] = xs;
          maj_xs[i + 1] = xs;
        }
      }
    }

    // Now we construct the energy grid and majorant
    EGrid = std::make_shared<pndl::EnergyGrid>(egrid, settings::ngroups);
    Emaj = std::make_shared<pndl::CrossSection>(maj_xs, EGrid, 0);
  }
}

std::vector<BankedParticle> ImplicitLeakageDeltaTracker::transport(
    std::vector<Particle>& bank, bool noise,
    std::vector<BankedParticle>* noise_bank, const NoiseMaker* noise_maker) {
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread local storage
    ThreadLocalScores thread_scores;

// Transport all particles in for thread
#ifdef _OPENMP
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
        warning(mssg.str(), __FILE__, __LINE__);
        p.kill();
      }
      // Only make helper if we aren't lost, to make sure that material isn't
      // a nullptr
      MaterialHelper mat(trkr.material(), p.E());
      if (settings::use_urr_ptables) mat.set_urr_rand_vals(p.rng);

      // auto bound = trkr.boundary();
      while (p.is_alive()) {
        bool had_collision = false;
        auto maj_indx = EGrid->get_lower_index(p.E());
        double Emajorant =
            Emaj->evaluate(p.E(), maj_indx) + mat.Ew(p.E(), noise);
        p.set_Esmp(Emajorant);  // Sampling XS saved for cancellation
        auto bound = trkr.boundary();

        double d_coll = 0.;
        if (bound.boundary_type == BoundaryType::Vacuum) {
          // Calculate probability of leaking and not leaking
          const double P_leak = std::exp(-Emajorant * bound.distance);
          const double P_no_leak = 1. - P_leak;

          // Get the weight which flys all the way to boundary and leaks
          const double wgt_leak = p.wgt() * P_leak;
          const double wgt2_leak = p.wgt2() * P_leak;

          // Get the weight which will have a collision
          const double wgt_collides = p.wgt() * P_no_leak;
          const double wgt2_collides = p.wgt2() * P_no_leak;

          // Score implicit leackage
          thread_scores.leakage_score += wgt_leak;

          // Score implicit migration area
          Position r_leak = p.r() + bound.distance * p.u();
          thread_scores.mig_score +=
              wgt_leak * (r_leak - p.r_birth()) * (r_leak - p.r_birth());

          // Sample flight distance
          d_coll = -std::log(1. - P_no_leak * RNG::rand(p.rng)) / Emajorant;

          // Score the TLE for the portion which would have leaked
          p.set_weight(wgt_leak);
          p.set_weight2(wgt2_leak);
          tallies->score_flight(p, bound.distance, mat, settings::converged);

          // Reduce weight
          p.set_weight(wgt_collides);
          p.set_weight2(wgt2_collides);

          // Score TLE for the portion which only goes to the collision site
          tallies->score_flight(p, d_coll, mat, settings::converged);
        } else {
          d_coll = RNG::exponential(p.rng, Emajorant);

          // Score track length tally for boundary distance.
          // This is here because flux-like tallies are allowed with DT.
          // No other quantity should be scored with a TLE, as an error
          // should have been thrown when building all tallies.
          tallies->score_flight(p, std::min(d_coll, bound.distance), mat,
                                settings::converged);
        }

        if (bound.distance < d_coll ||
            std::abs(bound.distance - d_coll) < 100. * SURFACE_COINCIDENT) {
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
              fatal_error(mssg.str(), __FILE__, __LINE__);
            }
          } else {
            fatal_error("Help me, how did I get here ?", __FILE__, __LINE__);
          }
        } else {
          // Update Position
          p.move(d_coll);
          trkr.move(d_coll);
          trkr.get_current();
          // Check if we are lost
          if (trkr.is_lost()) {
            std::stringstream mssg;
            mssg << "Particle " << p.history_id() << ".";
            mssg << p.secondary_id() << " has become lost.\n";
            mssg << "Previous valid coordinates: r = " << p.previous_r();
            mssg << ", u = " << p.previous_u() << ".\n";
            mssg << "Attempted to fly a distance of " << d_coll << " cm.\n";
            mssg << "Currently lost at r = " << trkr.r() << ", u = " << trkr.u()
                 << ".";
            fatal_error(mssg.str(), __FILE__, __LINE__);
          }
          mat.set_material(trkr.material(), p.E());

          // Get true cross section here
          double Et = mat.Et(p.E(), noise);

          if (Et - Emajorant > 1.E-10) {
            std::string mssg = "Total cross section excedeed majorant.";
            fatal_error(mssg, __FILE__, __LINE__);
          }

          if (RNG::rand(p.rng) < (Et / Emajorant)) {
            // Flag real collision
            had_collision = true;
          }
        }

        if (p.is_alive() && had_collision) {  // real collision
          collision(p, mat, thread_scores, noise, noise_maker);
          trkr.set_u(p.u());
          p.set_previous_collision_real();
        } else if (p.is_alive()) {  // Virtual collision
          p.set_previous_collision_virtual();
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
              fatal_error(mssg.str(), __FILE__, __LINE__);
            }
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
              warning(mssg, __FILE__, __LINE__);
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
  for (auto& p : bank) {
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
