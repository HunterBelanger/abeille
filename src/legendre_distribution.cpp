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
