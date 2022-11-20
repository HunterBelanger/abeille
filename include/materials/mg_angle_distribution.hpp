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
#ifndef MG_ANGLE_DISTRIBUTION_H
#define MG_ANGLE_DISTRIBUTION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <utils/rng.hpp>
#include <vector>

class MGAngleDistribution {
 public:
  // Isotropic Constructor
  MGAngleDistribution();

  // Anisotropic Constructor
  MGAngleDistribution(const std::vector<double>& mu,
                      const std::vector<double>& pdf,
                      const std::vector<double>& cdf);

  double sample_mu(pcg32& rng) const {
    const double xi = RNG::rand(rng);

    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    std::size_t l = std::distance(cdf_.begin(), cdf_it);
    if (xi == *cdf_it) return mu_[l];

    l--;

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    if (pdf_[l] == pdf_[l + 1]) return histogram_interp(xi, l);

    return linear_interp(xi, l);
  }

  double pdf(double mu) const {
    if (mu < min_value()) return pdf_.front();
    if (mu > max_value()) return pdf_.back();

    auto val_it = std::lower_bound(mu_.begin(), mu_.end(), mu);
    std::size_t l = std::distance(mu_.begin(), val_it);
    if (mu == *val_it) return pdf_[l];

    l--;

    const double m = (pdf_[l + 1] - pdf_[l]) / (mu_[l + 1] - mu_[l]);
    return m * (mu - mu_[l]) + pdf_[l];
  }

  double min_value() const { return mu_.front(); }

  double max_value() const { return mu_.back(); }

  const std::vector<double>& mu() const { return mu_; }

  const std::vector<double>& pdf() const { return pdf_; }

  const std::vector<double>& cdf() const { return cdf_; }

 private:
  std::vector<double> mu_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;

  double histogram_interp(double xi, std::size_t l) const {
    return mu_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp(double xi, std::size_t l) const {
    double m = (pdf_[l + 1] - pdf_[l]) / (mu_[l + 1] - mu_[l]);
    return mu_[l] +
           (1. / m) * (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) -
                       pdf_[l]);
  }
};

#endif
