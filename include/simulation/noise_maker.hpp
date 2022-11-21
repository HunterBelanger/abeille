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
#ifndef NOISE_MAKER_H
#define NOISE_MAKER_H

#include <yaml-cpp/yaml.h>

#include <materials/material_helper.hpp>
#include <memory>
#include <simulation/oscillation_noise_source.hpp>
#include <simulation/particle.hpp>
#include <simulation/vibration_noise_source.hpp>
#include <vector>

class NoiseMaker {
 public:
  NoiseMaker() {}

  void add_noise_source(const YAML::Node& snode);

  void add_noise_source(std::shared_ptr<VibrationNoiseSource> ns) {
    vibration_noise_sources_.push_back(ns);
  }

  void add_noise_source(std::shared_ptr<OscillationNoiseSource> ns) {
    oscillation_noise_sources_.push_back(ns);
  }

  std::size_t num_noise_sources() const {
    return vibration_noise_sources_.size() + oscillation_noise_sources_.size();
  }

  void sample_noise_source(Particle& p, MaterialHelper& mat, const double keff,
                           const double w) const;

 private:
  std::vector<std::shared_ptr<VibrationNoiseSource>> vibration_noise_sources_;
  std::vector<std::shared_ptr<OscillationNoiseSource>>
      oscillation_noise_sources_;

  bool is_inside(const Particle& p) const;
  std::complex<double> dEt(const Particle& p, double w) const;
  std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                          double w) const;
  std::unique_ptr<Material> make_fake_material(const Particle& p) const;

  void sample_noise_copy(Particle& p, MaterialHelper& mat,
                         const double w) const;

  void sample_vibration_noise_source(Particle& p, MaterialHelper& mat,
                                     const double keff, const double w) const;
  void sample_vibration_noise_fission(Particle& p, const Nuclide& nuclide,
                                      const MicroXSs& microxs,
                                      const std::complex<double>& dN_N,
                                      const double Etfake_Et, const double keff,
                                      const double w) const;
  void sample_vibration_noise_scatter(Particle& p, const Nuclide& nuclide,
                                      const MicroXSs& microxs,
                                      const std::complex<double>& dN_N,
                                      const double Etfake_Et,
                                      const double P_scatter) const;

  void sample_oscillation_noise_source(Particle& p, MaterialHelper& mat,
                                       const double keff, const double w) const;
  void sample_oscillation_noise_fission(Particle& p, const Nuclide& nuclide,
                                        const MicroXSs& microxs,
                                        const double keff,
                                        const double w) const;
  void sample_oscillation_noise_scatter(Particle& p, const Nuclide& nuclide,
                                        const MicroXSs& microxs,
                                        const double P_scatter,
                                        const double w) const;
};

#endif
