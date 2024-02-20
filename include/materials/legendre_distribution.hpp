/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef LEGENDRE_DISTRIBUTION_H
#define LEGENDRE_DISTRIBUTION_H

#include <materials/mg_angle_distribution.hpp>
#include <utils/rng.hpp>

#include <array>
#include <cmath>
#include <vector>

class LegendreDistribution {
 public:
  // Isotropic Constructor
  LegendreDistribution();
  // Anisotropic Constructor
  LegendreDistribution(const std::vector<double>& a);

  double sample_mu(RNG& rng) const {
    // Check for isotropic scattering
    if (a_.size() == 1) {
      return 2. * rng() - 1.;
    }

    double mu = -2.;
    while (true) {
      // First, sample which term to use, then sample the angle
      const double xi = rng();
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
      if (rng() < P_accept) break;
    }

    return mu;
  }

  double pdf(double mu) const {
    double p = 0;
    for (unsigned l = 0; l < a_.size(); l++) {
      p += a_[l] * legendre(l, mu);
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

  double sample_term_1(RNG& rng) const {
    const double xi = rng();
    return std::sin(PI * (xi - 0.5));
  }

  double sample_term_2(RNG& rng) const {
    const double xi = rng();
    const double a = 1.;
    const double b = 2.;
    const double c = 1. - (4. * xi);
    return (-b + std::sqrt(b * b - (4. * a * c))) / (2. * a);
  }

  double sample_term_3(RNG& rng) const {
    const double xi = rng();
    return 2. * xi - 1.;
  }

  double sample_term_4(RNG& rng) const {
    const double xi = rng();
    const double a = 1.;
    const double b = -2.;
    const double c = (4. * xi) - 3.;
    return (-b - std::sqrt(b * b - (4. * a * c))) / (2. * a);
  }

  // Copied from PNDL legendre.cpp
  double legendre(unsigned n, double x) const {
    switch (n) {
      case 0:
        return 1.;
        break;
      case 1:
        return x;
        break;
      case 2:
        return 0.5 * (3. * x * x - 1.);
        break;
      case 3:
        return 0.5 * (5. * x * x * x - 3. * x);
        break;
      case 4: {
        const double x_2 = x * x;
        return 0.125 * (35. * x_2 * x_2 - 30. * x_2 + 3.);
        break;
      }
      default: {
        // n >= 5, so start get getting p3 and p4
        const double x_2 = x * x;
        const double x_3 = x_2 * x;
        const double x_4 = x_3 * x;
        double p3 = 0.5 * (5. * x_3 - 3. * x);
        double p4 = 0.125 * (35. * x_4 - 30. * x_2 + 3.);

        // Iterate up to the desired Legendre order.
        unsigned l = 4;
        while (l < n) {
          std::swap(p3, p4);
          p4 = ((2. * l + 1.) * x * p3 - l * p4) / (l + 1.);
          l++;
        }
        return p4;

        break;
      }
    }
  }

  void initialize_values();
};

#endif
