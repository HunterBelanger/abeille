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
#include <materials/mg_angle_distribution.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <algorithm>
#include <sstream>

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
