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
#ifndef MG_ANGLE_DISTRIBUTION_H
#define MG_ANGLE_DISTRIBUTION_H

#include <utils/rng.hpp>

#include <algorithm>
#include <array>
#include <cmath>
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
    std::size_t l =
        static_cast<std::size_t>(std::distance(cdf_.begin(), cdf_it));
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
    std::size_t l =
        static_cast<std::size_t>(std::distance(mu_.begin(), val_it));
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
