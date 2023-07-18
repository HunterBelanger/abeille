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

#include <pcg_random.hpp>

#include <random>

class RNG {
 public:
  // -----------------------------------------------------------------------
  // rand()
  //   This advances the pcg32 generator and produces a double
  //   over the interval [0,1).
  // -----------------------------------------------------------------------
  static double rand(pcg32& rng) { return unit_dist(rng); }

  // -----------------------------------------------------------------------
  // uniform(double a, double b)
  //   Returns a random double with a uniform distribution over the
  //   interval [a, b)
  //
  //   f(x|a,b) = 1 / (b - a)
  // -----------------------------------------------------------------------
  static double uniform(pcg32& rng, double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // normal(double mu, double sigma)
  //   Returns a random double from the normal (Gaussian) distribution,
  //   defined by average value mu, and std sigma
  //
  //   f(x|mu,sigma) = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x - mu)/sigma)^2)
  // -----------------------------------------------------------------------
  static double normal(pcg32& rng, double mu, double sigma) {
    std::normal_distribution<double> dist(mu, sigma);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // exponential(double lambda)
  //   Returns a random double from the exponential distribution,
  //   defined by the average value 1/lambda.
  //
  //   f(x|lambda) = lambda * exp(-lambda * x)
  // -----------------------------------------------------------------------
  static double exponential(pcg32& rng, double lambda) {
    if (lambda == 0.) return INF;

    std::exponential_distribution<double> dist(lambda);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // discrete(std::vector<double> weights)
  //   Returns an integer over range [0,weights.size() - 1]
  //   where probability of each integer is defined by the weight
  //
  //   P(i|w_0, w_1, ... w_k) = w_i / Sum[j = 1 to k](w_j)
  // -----------------------------------------------------------------------
  static int discrete(pcg32& rng, const double* begin, const double* end) {
    std::discrete_distribution<int> dist(begin, end);
    return dist(rng);
  }

  static int discrete(pcg32& rng, const std::vector<double>& weights) {
    std::discrete_distribution<int> dist(weights.begin(), weights.end());
    return dist(rng);
  }

 private:
  static std::uniform_real_distribution<double> unit_dist;
};  // RNG

#endif
