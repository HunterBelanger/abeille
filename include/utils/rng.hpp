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

#include <cmath>
#include <span>
#include <vector>

class RNG {
 public:
  RNG() : engn() {}
  RNG(pcg32::state_type s) : engn(s) {}
  explicit RNG(pcg32 s) : engn(s) {}

  void seed(std::uint64_t s) { engn.seed(s); }

  void advance(std::uint64_t n) { engn.advance(n); }

  void set_stream(std::uint32_t s) { engn.set_stream(s); }

  double rand() { return pcg_to_double(engn()); }

  double operator()() { return this->rand(); }

  double uniform(double a, double b) { return (b - a) * this->rand() + a; }

  double exponential(double lambda) {
    if (lambda < 0.) {
      fatal_error("Lambda must be >= 0.");
    }

    if (lambda == 0.) return INF;

    return -std::log(rand()) / lambda;
  }

  std::size_t discrete(std::span<const double> weights) {
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

    // If sum of weights is zero, we can't sample a value. The STL
    // distributions seems to just return 0 in this case.
    if (norm == 0.) {
      return 0;
    }

    // Sample the index
    const double xi = rand() * norm;
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

  std::size_t discrete(const std::vector<double>& weights) {
    std::span<const double> spn(weights.begin(), weights.end());
    return this->discrete(spn);
  }

  std::uint64_t operator-(const RNG& rhs) const { return engn - rhs.engn; }

  static double max() { return pcg_to_double(pcg32::max()); }

  static double min() { return pcg_to_double(pcg32::min()); }

 private:
  pcg32 engn;

  static double pcg_to_double(pcg32::result_type val) {
    return std::ldexp(static_cast<double>(val), -32);
  }
};

#endif
