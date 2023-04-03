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
#include <cmath>
#include <memory>
#include <simulation/cone.hpp>
#include <utils/error.hpp>

Cone::Cone(Direction u, double aperture)
    : direction_(u), aperture_(std::cos(aperture)) {}

Direction Cone::sample(pcg32& rng) const {
  // Sample mu in [aperture_, 1]
  double mu = (1. - aperture_) * RNG::rand(rng) + aperture_;

  // Sample phi in [0, 2pi]
  double phi = 2. * PI * RNG::rand(rng);

  return rotate_direction(direction_, mu, phi);
}

double Cone::aperture() const { return std::acos(aperture_); }

std::shared_ptr<Cone> make_cone_distribution(const YAML::Node& node) {
  if (!node["direction"] || !node["direction"].IsSequence() ||
      !(node["direction"].size() == 3)) {
    fatal_error("No valid direction entry for cone distribution.");
  }

  double x = node["direction"][0].as<double>();
  double y = node["direction"][1].as<double>();
  double z = node["direction"][2].as<double>();

  Direction u(x, y, z);

  if (!node["aperture"] || !node["aperture"].IsScalar()) {
    fatal_error("No valid aperture entry for cone distribution.");
  }
  double aperture = node["aperture"].as<double>();

  return std::make_shared<Cone>(u, aperture);
}
