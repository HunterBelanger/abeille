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
#ifndef VIBRATION_NOISE_SOURCE_H
#define VIBRATION_NOISE_SOURCE_H

#include <complex>
#include <memory>
#include <simulation/noise_source.hpp>
#include <unordered_map>
#include <vector>

// Pure virtual interface to allow for the sampling of neutron noise particles.
class VibrationNoiseSource : public NoiseSource {
 public:
  struct NuclideInfo {
    uint32_t id;
    double concentration;
  };

  VibrationNoiseSource() : nuclides_(), nuclide_info_() {}
  virtual ~VibrationNoiseSource() = default;

  virtual std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                                  double w) const = 0;

  const std::vector<uint32_t>& nuclides() const { return nuclides_; }

  const std::unordered_map<uint32_t, NuclideInfo>& nuclide_info() const {
    return nuclide_info_;
  }

 protected:
  std::vector<uint32_t> nuclides_;
  std::unordered_map<uint32_t, NuclideInfo> nuclide_info_;
};

#endif
