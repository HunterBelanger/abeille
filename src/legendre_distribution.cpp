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
#include <materials/legendre_distribution.hpp>
#include <utils/constants.hpp>
#include <vector>

LegendreDistribution::LegendreDistribution()
    : a_({0.5}),
      cdf_({0., 0., 0., 0.}),
      C1_(0.),
      C2_(0.),
      C3_(0.),
      C4_(0.),
      H_(0.) {
  initialize_values();
}

LegendreDistribution::LegendreDistribution(const std::vector<double>& a)
    : a_(a),
      cdf_({0., 0., 0., 0.}),
      C1_(0.),
      C2_(0.),
      C3_(0.),
      C4_(0.),
      H_(0.) {
  // Add the 0th moment, which is always 1
  a_.insert(a_.begin(), 1.);

  // Now we need to turn all of the legendre moments into the
  // coefficients required by the algorithm in Lux and Koblinger.
  // To do this, we multiply all elements a_[i] by (2i + 1)/2.
  for (std::size_t l = 0; l < a_.size(); l++) {
    a_[l] *= ((2. * static_cast<double>(l) + 1.) / 2.);
  }

  initialize_values();
}

void LegendreDistribution::initialize_values() {
  // Calculate C1
  for (std::size_t l = 2; l < a_.size(); l++) {
    C1_ += std::abs(a_[l]) * std::sqrt(static_cast<double>(l));
  }
  C1_ *= std::sqrt(2. / PI);

  // Calculate C2
  if (a_.size() > 1) C2_ = std::max(a_[1], 0.);

  // Calculate C3
  if (a_.size() > 1) {
    if (a_[0] >= a_[1] && a_[1] > 0.)
      C3_ = a_[0] - a_[1];
    else if (-a_[0] <= a_[1] && a_[1] < 0.)
      C3_ = a_[0] + a_[1];
  }

  // Calculate C4
  if (a_.size() > 1) C4_ = std::max(-a_[1], 0.);

  // Calculate H
  H_ = PI * C1_ + 2. * (C2_ + C3_ + C4_);

  // Calculate coefficients
  std::array<double, 4> probs;
  probs.fill(0.);
  probs[0] = PI * C1_ / H_;
  probs[1] = 2. * C2_ / H_;
  probs[2] = 2. * C3_ / H_;
  probs[3] = 2. * C4_ / H_;

  cdf_[0] = probs[0];
  cdf_[1] = cdf_[0] + probs[1];
  cdf_[2] = cdf_[1] + probs[2];
  cdf_[3] = cdf_[2] + probs[3];

  // Make sure cdf_ is normalized
  for (std::size_t i = 0; i < cdf_.size(); i++) {
    cdf_[i] /= cdf_[3];
  }
}

MGAngleDistribution LegendreDistribution::linearize() const {
  std::vector<double> mu{-1., 1.};
  std::vector<double> p;
  p.push_back(pdf(-1.));
  p.push_back(pdf(1.));

  // Bisect intervals until we are linearly interpolable.
  std::size_t i = 0;
  while (i < (mu.size() - 1)) {
    // Get the mid-point value
    double mu_mid = 0.5 * (mu[i] + mu[i + 1]);

    // Get interpolated and real pdf
    double p_interp = 0.5 * (p[i] + p[i + 1]);
    double p_real = pdf(mu_mid);

    // Check tolerance
    double rel_diff = std::abs(p_interp - p_real) / p_real;
    if (rel_diff > TOLERANCE) {
      // We need to add a new point
      auto ip =
          mu.begin() + static_cast<std::vector<double>::difference_type>(i) + 1;
      auto pp =
          p.begin() + static_cast<std::vector<double>::difference_type>(i) + 1;
      mu.insert(ip, mu_mid);
      p.insert(pp, p_real);
    } else {
      i++;
    }
  }

  // Initialize all zero CDF vector.
  std::vector<double> cdf(mu.size(), 0.);

  // Trapezoid rule for integral of PDF, stored in the CDF.
  for (std::size_t i = 0; i < mu.size() - 1; i++) {
    cdf[i + 1] = ((mu[i + 1] - mu[i]) * 0.5 * (p[i + 1] + p[i])) + cdf[i];
  }

  // Normalize the PDF and CDF by the integral of last value in the
  // CDF, which should be integral(PDF(x), -1, 1), and this then
  // ensures we have proper normalization.
  const double norm = cdf.back();
  for (std::size_t i = 0; i < cdf.size(); i++) {
    p[i] /= norm;
    cdf[i] /= norm;
  }

  // Construct and return the linearly interpolable distribution.
  return MGAngleDistribution(mu, p, cdf);
}
