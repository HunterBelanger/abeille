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
#ifndef LEGENDRE_DISTRIBUTION_H
#define LEGENDRE_DISTRIBUTION_H
#include <array>
#include <cmath>
#include <materials/mg_angle_distribution.hpp>
#include <utils/rng.hpp>
#include <vector>

class LegendreDistribution {
 public:
  // Isotropic Constructor
  LegendreDistribution();
  // Anisotropic Constructor
  LegendreDistribution(const std::vector<double>& a);

  double sample_mu(pcg32& rng) const {
    // Check for isotropic scattering
    if (a_.size() == 1) {
      return 2. * RNG::rand(rng) - 1.;
    }

    double mu = -2.;
    while (true) {
      // First, sample which term to use, then sample the angle
      const double xi = RNG::rand(rng);
      if (xi <= cdf_[0]) {
        mu = sample_term_1(rng);
      } else if (xi <= cdf_[1]) {
        mu = sample_term_2(rng);
      } else if (xi <= cdf_[2]) {
        mu = sample_term_3(rng);
      } else {
        mu = sample_term_4(rng);
      }

      // Make sure mu is within [-1, 1]
      clamp_mu(mu);

      // Check if we accept or reject
      const double P_accept = pdf(mu) / h(mu);
      if (RNG::rand(rng) < P_accept) break;
    }

    return mu;
  }

  double pdf(double mu) const {
    double p = 0;
    for (std::size_t l = 0; l < a_.size(); l++) {
      p += a_[l] * std::legendre(l, mu);
    }
    return p;
  }

  double h(double mu) const {
    double h = 0.;
    h += C1_ / std::sqrt(1. - (mu * mu));
    h += C2_ * (1. + mu);
    h += C3_;
    h += C4_ * (1. - mu);
    return h;
  }

  const std::vector<double>& a() const { return a_; }

  void set_moment(std::size_t l, double coeff) {
    // No touching the 0th moment. Then we mess up normalization.
    if (l == 0) return;

    // Add moments which have coeffs of 0 to get
    // enough entries if needed.
    a_.resize(l + 1, 0.);

    a_[l] = coeff * (2. * static_cast<double>(l) + 1.) / 2.;

    initialize_values();
  }

  bool positive_over_domain() const {
    double mu = -1.;
    constexpr double d_mu = 0.01;

    for (std::size_t i = 0; i <= 200; i++) {
      if (pdf(mu) < 0.) {
        return false;
      }
      mu += d_mu;
    }

    return true;
  }

  MGAngleDistribution linearize() const;

 private:
  std::vector<double> a_;
  std::array<double, 4> cdf_;
  double C1_, C2_, C3_, C4_, H_;

  void clamp_mu(double& mu) const {
    if (mu > 1.)
      mu = 1.;
    else if (mu < -1.)
      mu = -1.;
  }

  double sample_term_1(pcg32& rng) const {
    const double xi = RNG::rand(rng);
    return std::sin(PI * (xi - 0.5));
  }

  double sample_term_2(pcg32& rng) const {
    const double xi = RNG::rand(rng);
    const double a = 1.;
    const double b = 2.;
    const double c = 1. - (4. * xi);
    return (-b + std::sqrt(b * b - (4. * a * c))) / (2. * a);
  }

  double sample_term_3(pcg32& rng) const {
    const double xi = RNG::rand(rng);
    return 2. * xi - 1.;
  }

  double sample_term_4(pcg32& rng) const {
    const double xi = RNG::rand(rng);
    const double a = 1.;
    const double b = -2.;
    const double c = (4. * xi) - 3.;
    return (-b - std::sqrt(b * b - (4. * a * c))) / (2. * a);
  }

  void initialize_values();
};

#endif
