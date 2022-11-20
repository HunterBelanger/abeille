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
#ifndef FLAT_VIBRATION_NOISE_SOURCE_H
#define FLAT_VIBRATION_NOISE_SOURCE_H

#include <complex>
#include <materials/material.hpp>
#include <memory>
#include <simulation/vibration_noise_source.hpp>
#include <unordered_map>

class FlatVibrationNoiseSource : public VibrationNoiseSource {
 public:
  enum class Basis { X, Y, Z };

 public:
  FlatVibrationNoiseSource(Position low, Position hi, Basis basis,
                           std::shared_ptr<Material> mat_pos,
                           std::shared_ptr<Material> mat_neg,
                           double angular_frequency);

  bool is_inside(const Position& r) const override final;

  std::complex<double> dEt(const Position& r, double E,
                           double w) const override final;

  std::complex<double> dEt_Et(const Position& r, double E,
                              double w) const override final;

  std::complex<double> dN(const Position& r, uint32_t nuclide_id,
                          double w) const override final;

 private:
  Position low_, hi_;
  Basis basis_;
  std::shared_ptr<Material> material_pos_;  // Material on the positive side
  std::shared_ptr<Material> material_neg_;  // Material on the negative side
  double x0_;                               // Midpoint of interface
  double w0_;                               // Angular frequency of vibration
  double eps_;                              // Magnitude of oscillation
  std::unordered_map<uint32_t, double> Delta_N;  // Contians N_neg - N_pos

  // Function to get xs
  double Et(double x, double E) const;
  // Function to get Ei_neg(E) - Ei_pos(E)
  double Delta_Et(double E) const;

  std::complex<double> C_R(uint32_t n, double x) const;
  std::complex<double> C_L(uint32_t n, double x) const;

  double get_x(const Position& r) const;

  bool negative_material(double x) const;

  double get_nuclide_concentration(const Material& mat,
                                   uint32_t nuclide_id) const;

  static constexpr std::complex<double> i{0., 1.};
};

std::shared_ptr<VibrationNoiseSource> make_flat_vibration_noise_source(
    const YAML::Node& snode);

#endif
