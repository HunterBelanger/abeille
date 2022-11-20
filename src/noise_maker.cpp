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
#include <algorithm>
#include <complex>
#include <iterator>
#include <materials/material.hpp>
#include <materials/nuclide.hpp>
#include <memory>
#include <simulation/flat_vibration_noise_source.hpp>
#include <simulation/noise_maker.hpp>
#include <simulation/square_oscillation_noise_source.hpp>
#include <utils/error.hpp>
#include <utils/rng.hpp>
#include <vector>

void NoiseMaker::add_noise_source(const YAML::Node& snode) {
  // Make sure it is a map
  if (!snode.IsMap()) {
    std::string mssg = "Noise source entry must be a map.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get the type of the source
  if (!snode["type"] || !snode["type"].IsScalar()) {
    std::string mssg = "No valid type provided to noise source entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::string type = snode["type"].as<std::string>();

  if (type == "flat-vibration") {
    this->add_noise_source(make_flat_vibration_noise_source(snode));
  } else if (type == "square-oscillation") {
    this->add_noise_source(make_square_oscillation_noise_source(snode));
  } else {
    std::string mssg = "Invalid noise source type " + type + ".";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

std::complex<double> NoiseMaker::dEt(const Particle& p, double w) const {
  std::complex<double> dEt_to_return{0., 0.};

  for (const auto& ns : vibration_noise_sources_) {
    // Only contribute if we are inside the source region.
    if (ns->is_inside(p.r())) {
      dEt_to_return += ns->dEt(p.r(), p.E(), w);
    }
  }

  for (const auto& ns : oscillation_noise_sources_) {
    // Only contribute if we are inside the source region.
    if (ns->is_inside(p.r())) {
      dEt_to_return += ns->dEt(p.r(), p.E(), w);
    }
  }

  return dEt_to_return;
}

std::complex<double> NoiseMaker::dN(const Position& r, uint32_t nuclide_id,
                                    double w) const {
  std::complex<double> dN_to_return{0., 0.};

  for (const auto& ns : vibration_noise_sources_) {
    // Only contribute if we are inside the source region.
    if (ns->is_inside(r)) {
      dN_to_return += ns->dN(r, nuclide_id, w);
    }
  }

  return dN_to_return;
}

bool NoiseMaker::is_inside(const Particle& p) const {
  // To be inside, it sufices to be inside only one noise source.
  for (const auto& ns : vibration_noise_sources_) {
    if (ns->is_inside(p.r())) return true;
  }

  for (const auto& ns : oscillation_noise_sources_) {
    if (ns->is_inside(p.r())) return true;
  }

  return false;
}

std::unique_ptr<Material> NoiseMaker::make_fake_material(
    const Particle& p) const {
  std::vector<uint32_t> nuclides_list;
  std::vector<std::shared_ptr<VibrationNoiseSource>> sources_list;

  // Go through all noise sources
  for (const auto& ns : vibration_noise_sources_) {
    // Check if we are inside it
    if (ns->is_inside(p.r())) {
      // Save noise source
      sources_list.push_back(ns);

      // Add these nuclides to the nuclides list
      std::vector<uint32_t> tmp;
      std::set_union(nuclides_list.begin(), nuclides_list.end(),
                     ns->nuclides().begin(), ns->nuclides().end(),
                     std::back_inserter(tmp));
      nuclides_list = tmp;
    }
  }

  // nuclides_list should now contain all of the nuclides which are
  // involved in producing noise. We need to calculate the concentration
  // for each of them in the fake material now.
  std::vector<double> concentrations(nuclides_list.size(), 0.);

  const double num_sources = static_cast<double>(sources_list.size());
  for (std::size_t i = 0; i < nuclides_list.size(); i++) {
    const uint32_t nuclide_id = nuclides_list[i];
    double conc_sum = 0.;
    for (const auto& ns : sources_list) {
      conc_sum += ns->nuclide_info().at(nuclide_id).concentration;
    }
    concentrations[i] = conc_sum / num_sources;
  }

  // Now we make a new material, and add all of the components.
  std::unique_ptr<Material> fake_mat = std::make_unique<Material>();

  for (std::size_t i = 0; i < nuclides_list.size(); i++) {
    const uint32_t& nuclide_id = nuclides_list[i];
    auto nuclide = nuclides.at(nuclide_id);
    fake_mat->add_component(concentrations[i], nuclide);
  }

  return fake_mat;
}

void NoiseMaker::sample_noise_copy(Particle& p, MaterialHelper& mat,
                                   const double w) const {
  const std::complex<double> dEt_Et = this->dEt(p, w) / mat.Et(p.E());

  std::complex<double> weight_copy{p.wgt(), p.wgt2()};

  // Negative needed for source sampling
  weight_copy *= -dEt_Et;

  // Create and return the neutron
  BankedParticle p_noise{p.r(),
                         p.u(),
                         p.E(),
                         weight_copy.real(),
                         weight_copy.imag(),
                         p.history_id(),
                         p.daughter_counter(),
                         p.previous_collision_virtual(),
                         p.previous_r(),
                         p.previous_u(),
                         p.previous_E(),
                         p.E(),
                         p.Esmp()};

  p.add_noise_particle(p_noise);
}

void NoiseMaker::sample_vibration_noise_fission(
    Particle& p, const std::shared_ptr<Nuclide>& nuclide,
    const MicroXSs& microxs, const std::complex<double>& dN_N,
    const double Etfake_Et, const double keff, const double w) const {
  // First, check if the nuclide is fissile. If it isn't, we
  // just return.
  if (nuclide->fissile() == false) return;

  // Now we get the number of noise neutrons to produce.
  const double k_abs = microxs.nu_total * microxs.fission / microxs.total;
  const int n_new =
      static_cast<int>(std::floor(k_abs / keff + RNG::rand(p.rng)));

  // Calculate the probability of a delayed neutron
  const double P_delayed = microxs.nu_delayed / microxs.nu_total;

  for (int n = 0; n < n_new; n++) {
    auto finfo = nuclide->sample_fission(p.E(), p.u(), microxs.energy_index,
                                         P_delayed, p.rng);

    // Sample a banked noise particle from fission.
    BankedParticle bnp{p.r(),
                       finfo.direction,
                       finfo.energy,
                       p.wgt(),
                       p.wgt2(),
                       p.history_id(),
                       p.daughter_counter(),
                       p.previous_collision_virtual(),
                       p.previous_r(),
                       p.previous_u(),
                       p.previous_E(),
                       p.E(),
                       p.Esmp()};

    if (finfo.delayed) {
      std::complex<double> wgt_cmpx{bnp.wgt, bnp.wgt2};
      double lambda = finfo.precursor_decay_constant;
      double denom = (lambda * lambda) + (w * w);
      std::complex<double> mult{lambda * lambda / denom, -lambda * w / denom};
      wgt_cmpx *= mult;
      bnp.wgt = wgt_cmpx.real();
      bnp.wgt2 = wgt_cmpx.imag();
    }

    // Apply the weight corrections
    std::complex<double> bnp_wgt{bnp.wgt, bnp.wgt2};
    bnp_wgt *= dN_N;
    bnp_wgt *= Etfake_Et;
    bnp.wgt = bnp_wgt.real();
    bnp.wgt2 = bnp_wgt.imag();

    // Save the noise particle
    p.add_noise_particle(bnp);
  }
}

void NoiseMaker::sample_vibration_noise_scatter(
    Particle& p, const std::shared_ptr<Nuclide>& nuclide,
    const MicroXSs& microxs, const std::complex<double>& dN_N,
    const double Etfake_Et, const double P_scatter) const {
  // First, sample the scatter info from the nuclide
  ScatterInfo sinfo =
      nuclide->sample_scatter(p.E(), p.u(), microxs.energy_index, p.rng);

  // Make the noise particle without weights
  BankedParticle p_noise{p.r(),
                         sinfo.direction,
                         sinfo.energy,
                         0.,  // wgt
                         0.,  // wgt2
                         p.history_id(),
                         p.daughter_counter(),
                         p.previous_collision_virtual(),
                         p.previous_r(),
                         p.previous_u(),
                         p.previous_E(),
                         p.E(),
                         p.Esmp()};

  std::complex<double> wgt{p.wgt(), p.wgt2()};
  wgt *= sinfo.yield;       // Multiply by scattering yield
  wgt *= P_scatter;         // Implicit capture
  wgt *= dN_N * Etfake_Et;  // Weight corrections for noise

  // Set the weight
  p_noise.wgt = wgt.real();
  p_noise.wgt2 = wgt.imag();

  // Save noise particle
  p.add_noise_particle(p_noise);
}

void NoiseMaker::sample_noise_source(Particle& p, MaterialHelper& mat,
                                     const double keff, const double w) const {
  // First, check if we are inside any noise source. If no, we can return
  // without making any noise particles.
  if (this->is_inside(p) == false) return;

  // First, we can go ahead and make the copy, as it's easiest to do.
  this->sample_noise_copy(p, mat, w);

  // Now from this point on, we sample the oscillation and vibration
  // parts separately
  this->sample_vibration_noise_source(p, mat, keff, w);

  this->sample_oscillation_noise_source(p, mat, keff, w);
}

void NoiseMaker::sample_oscillation_noise_source(Particle& p,
                                                 MaterialHelper& mat,
                                                 const double keff,
                                                 const double w) const {
  // First, check if we are inside any oscillation noise source.
  // If no, we can return without making any noise particles.
  bool inside_oscillation_source = false;
  for (const auto& ns : oscillation_noise_sources_) {
    if (ns->is_inside(p.r())) {
      inside_oscillation_source = true;
      break;
    }
  }
  if (inside_oscillation_source == false) {
    return;
  }

  // We now sample a nuclide for sampling the source particles.
  auto nuclide_data = mat.sample_nuclide(p.E(), p.rng);
  const MicroXSs& microxs = nuclide_data.second;
  const std::shared_ptr<Nuclide>& nuclide = nuclide_data.first;

  this->sample_oscillation_noise_fission(p, nuclide, microxs, keff, w);

  // Now we calculate the scatter probability to modify the weight
  // of the scatter noise source. This is the "implicit capture" correciton,
  // needed because we forced the sampling of the fission noise source.
  const double P_scatter = 1. - (microxs.absorption / microxs.total);

  this->sample_oscillation_noise_scatter(p, nuclide, microxs, P_scatter, w);
}

void NoiseMaker::sample_oscillation_noise_scatter(
    Particle& p, const std::shared_ptr<Nuclide>& nuclide,
    const MicroXSs& microxs, const double P_scatter, const double w) const {
  // First, sample the scatter info from the nuclide
  ScatterInfo sinfo =
      nuclide->sample_scatter(p.E(), p.u(), microxs.energy_index, p.rng);

  // Make the noise particle without weights
  BankedParticle p_noise{p.r(),
                         sinfo.direction,
                         sinfo.energy,
                         0.,  // wgt
                         0.,  // wgt2
                         p.history_id(),
                         p.daughter_counter(),
                         p.previous_collision_virtual(),
                         p.previous_r(),
                         p.previous_u(),
                         p.previous_E(),
                         p.E(),
                         p.Esmp()};

  std::complex<double> wgt{p.wgt(), p.wgt2()};
  wgt *= sinfo.yield;  // Multiply by scattering yield
  wgt *= P_scatter;    // Implicit capture

  // Get the complex weight factor
  std::complex<double> dE_E{0., 0.};
  if (sinfo.mt == 2) {
    // Go through and get all elastic scatter stuff
    for (const auto& ns : oscillation_noise_sources_) {
      if (ns->is_inside(p.r())) {
        dE_E += ns->dEelastic_Eelastic(p.r(), p.E(), w);
      }
    }
  } else {
    // Go through and get MT info
    // ONLY POSSIBLE IN CE MODE
    for (const auto& ns : oscillation_noise_sources_) {
      if (ns->is_inside(p.r())) {
        dE_E += ns->dEmt_Emt(sinfo.mt, p.r(), p.E(), w);
      }
    }
  }

  wgt *= dE_E;

  // Set the weight
  p_noise.wgt = wgt.real();
  p_noise.wgt2 = wgt.imag();

  // Save noise particle
  p.add_noise_particle(p_noise);
}

void NoiseMaker::sample_oscillation_noise_fission(
    Particle& p, const std::shared_ptr<Nuclide> nuclide,
    const MicroXSs& microxs, const double keff, const double w) const {
  // First, check if the nuclide is fissile. If it isn't, we
  // just return.
  if (nuclide->fissile() == false) return;

  // Now we get the number of noise neutrons to produce.
  const double k_abs = microxs.nu_total * microxs.fission / microxs.total;
  const int n_new =
      static_cast<int>(std::floor(k_abs / keff + RNG::rand(p.rng)));

  // Calculate the probability of a delayed neutron
  const double P_delayed = microxs.nu_delayed / microxs.nu_total;

  // Get the weight modifier
  std::complex<double> dEf_Ef = {0., 0.};
  for (const auto& ns : oscillation_noise_sources_) {
    if (ns->is_inside(p.r())) {
      dEf_Ef += ns->dEf_Ef(p.r(), p.E(), w);
    }
  }

  for (int i = 0; i < n_new; i++) {
    auto finfo = nuclide->sample_fission(p.E(), p.u(), microxs.energy_index,
                                         P_delayed, p.rng);

    BankedParticle bnp{p.r(),
                       finfo.direction,
                       finfo.energy,
                       p.wgt(),
                       p.wgt2(),
                       p.history_id(),
                       p.daughter_counter(),
                       p.previous_collision_virtual(),
                       p.previous_r(),
                       p.previous_u(),
                       p.previous_E(),
                       p.E(),
                       p.Esmp()};

    if (finfo.delayed) {
      std::complex<double> wgt_cmpx{bnp.wgt, bnp.wgt2};
      double lambda = finfo.precursor_decay_constant;
      double denom = (lambda * lambda) + (w * w);
      std::complex<double> mult{lambda * lambda / denom, -lambda * w / denom};
      wgt_cmpx *= mult;
      bnp.wgt = wgt_cmpx.real();
      bnp.wgt2 = wgt_cmpx.imag();
    }

    // Modify weight
    std::complex<double> fnp_weight{bnp.wgt, bnp.wgt2};
    fnp_weight *= dEf_Ef;
    bnp.wgt = fnp_weight.real();
    bnp.wgt2 = fnp_weight.imag();

    // Save BankedParticle
    p.add_noise_particle(bnp);
  }
}

void NoiseMaker::sample_vibration_noise_source(Particle& p, MaterialHelper& mat,
                                               const double keff,
                                               const double w) const {
  // First, check if we are inside any vibration noise source.
  // If no, we can return without making any noise particles.
  bool inside_vibration_source = false;
  for (const auto& ns : vibration_noise_sources_) {
    if (ns->is_inside(p.r())) {
      inside_vibration_source = true;
      break;
    }
  }
  if (inside_vibration_source == false) {
    return;
  }

  // Now, we need to get our ficticious material, which represents a
  // sort of homogenization of all materials involved in all of the
  // noise regions where the particle is currently located.
  std::unique_ptr<Material> fake_material = make_fake_material(p);
  MaterialHelper fake_mat(fake_material.get(), p.E());

  // We now sample a nuclide from this fake material.
  auto nuclide_data = fake_mat.sample_nuclide(p.E(), p.rng);
  const MicroXSs& microxs = nuclide_data.second;
  const std::shared_ptr<Nuclide>& nuclide = nuclide_data.first;
  // Id of sampled nuclide
  const uint32_t nuclide_id = nuclide->id();
  // Concentration of sampled nuclide in fake material
  const double N = microxs.concentration;

  // Ratio of dN / N, to be used as part of the weight modifier.
  const std::complex<double> dN_N = this->dN(p.r(), nuclide_id, w) / N;

  // Ratio of Et in the fake material, to Et in the actual material.
  // This will also be used as a weight modifier.
  const double Etfake_Et = fake_mat.Et(p.E()) / mat.Et(p.E());

  // Sample the fission noise source particles
  this->sample_vibration_noise_fission(p, nuclide, microxs, dN_N, Etfake_Et,
                                       keff, w);

  // Now we calculate the scatter probability to modify the weight
  // of the scatter noise source. This is the "implicit capture" correciton,
  // needed because we forced the sampling of the fission noise source.
  const double P_scatter = 1. - (microxs.absorption / microxs.total);

  // Sample the scatter noise source particles
  this->sample_vibration_noise_scatter(p, nuclide, microxs, dN_N, Etfake_Et,
                                       P_scatter);
}
