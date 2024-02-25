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
#ifndef RNG_H
#define RNG_H

#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <pcg_random.hpp>
#include <gsl/gsl-lite.hpp>

#include <cmath>
#include <vector>

class RNG {
 public:
  // -----------------------------------------------------------------------
  // rand()
  //   This advances the pcg32 generator and produces a double
  //   over the interval [0,1).
  // -----------------------------------------------------------------------
  static double rand(pcg32& rng) {
    return pcg_to_double(rng());
  }

  // -----------------------------------------------------------------------
  // uniform(double a, double b)
  //   Returns a random double with a uniform distribution over the
  //   interval [a, b)
  //
  //   f(x|a,b) = 1 / (b - a)
  // -----------------------------------------------------------------------
  static double uniform(pcg32& rng, double a, double b) {
    return (b - a) * rand(rng) + a;
  }

  // -----------------------------------------------------------------------
  // exponential(double lambda)
  //   Returns a random double from the exponential distribution,
  //   defined by the average value 1/lambda.
  //
  //   f(x|lambda) = lambda * exp(-lambda * x)
  // -----------------------------------------------------------------------
  static double exponential(pcg32& rng, double lambda) {
    if (lambda < 0.) {
      fatal_error("Lambda must be >= 0.");
    }

    if (lambda == 0.) return INF;

    return -std::log(rand(rng)) / lambda;
  }

  // -----------------------------------------------------------------------
  // discrete(std::vector<double> weights)
  //   Returns an integer over range [0,weights.size() - 1]
  //   where probability of each integer is defined by the weight
  //
  //   P(i|w_0, w_1, ... w_k) = w_i / Sum[j = 1 to k](w_j)
  // -----------------------------------------------------------------------
  static std::size_t discrete(pcg32& rng, gsl::span<const double> weights) {
    // If weights is empty, we can't sample an index. In this case, we
    // return 0 like for the STL's discrete_distribution.
    if (weights.empty()) {
      return 0;
    }

    // Calculate the sum of all weights
    double norm = 0.;
    for (const auto& val : weights) {
      // If a weight is < 0, the distribution is not well defined.
      if (val < 0.) {
        fatal_error("Weights must be >= 0.");
      }

      norm += val;
    }

    // If sum of weights is zero, we can't sample a value
    if (norm == 0.) {
      return 0;
      //fatal_error("Sum of all weights must be > 0.");
    }

    // Sample the index
    const double xi = rand(rng) * norm;
    double sum = 0.;
    for (std::size_t i = 0; i < weights.size(); i++) {
      sum += weights[i];

      if (xi < sum) {
        return i;
      }
    }

    // Should never get here, but return last index
    return weights.size() - 1;
  }

  static std::size_t discrete(pcg32& rng, const double* begin, const double* end) {
    return discrete(rng, gsl::span<const double>(begin, end));
  }

  static std::size_t discrete(pcg32& rng, const std::vector<double>& weights) {
    gsl::span<const double> spn(&weights[0], &weights[0]+weights.size());
    return discrete(rng, spn);
  }

  static double pcg_to_double(pcg32::result_type val) {
    return std::ldexp(static_cast<double>(val), -32);
  }
};  // RNG

#endif
