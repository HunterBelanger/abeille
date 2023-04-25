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
#include <materials/mg_angle_distribution.hpp>
#include <sstream>
#include <utils/constants.hpp>
#include <utils/error.hpp>

MGAngleDistribution::MGAngleDistribution()
    : mu_({-1., 1.}), pdf_({0.5, 0.5}), cdf_({0., 1.}) {}

MGAngleDistribution::MGAngleDistribution(const std::vector<double>& mu,
                                         const std::vector<double>& pdf,
                                         const std::vector<double>& cdf)
    : mu_(mu), pdf_(pdf), cdf_(cdf) {
  // Make sure good mu bounds
  if (mu_.front() < -1.) {
    fatal_error("Angle limit less than -1.");
  }

  if (mu_.front() > 1.) {
    fatal_error("Angle limit greater than 1.");
  }

  // Make sure mu is sorted
  if (std::is_sorted(mu_.begin(), mu_.end()) == false) {
    fatal_error("Mu values are not sorted.");
  }

  // Make sure PDF is positive
  for (const auto& p : pdf_) {
    if (p < 0.) {
      fatal_error("PDF is less than 0.");
    }
  }

  // Make sure CDF is positive
  for (const auto& c : cdf_) {
    if (c < 0.) {
      fatal_error("CDF is less than 0.");
    }
  }

  // Make sure CDF is sorted
  if (std::is_sorted(cdf_.begin(), cdf_.end()) == false) {
    fatal_error("CDF is not sorted.");
  }

  // Make sure CDF starts at 0, and ends at 1
  if (cdf_.front() != 0.) {
    fatal_error("First CDF value is not 0.");
  }

  if (cdf_.back() != 1.) {
    fatal_error("Last CDF value is not 1.");
  }
}
