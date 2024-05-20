#ifndef ZERNIKE_POLYNOMIAL_H
#define ZERNIKE_POLYNOMIAL_H

#include <utils/error.hpp>

#include <cmath>
#include <iostream>
#include <vector>

//===================================
// Zernike Polynomials from order 0 up to order.
//
class ZernikePolynomials {
 public:
  // enum class for the even and odd zernike
  enum class ZernikeType { Even, Odd };

  ZernikePolynomials(std::size_t order);

  double evaluate_zernike_at_order(const double x, const double theta,
                                   const std::size_t order) const;

  std::vector<double> evaluate_zernikes(const double r,
                                        const double theta) const;

 private:
  std::size_t order_, max_n_;
  std::vector<std::vector<double>> Zr_coefficients_;
  std::vector<std::size_t> m_;
  std::vector<std::size_t> n_;
  std::vector<ZernikeType> Zr_types_;

  // function for the factorial
  std::size_t factorial(std::size_t N) {
    if (N == 1 || N == 0)
      return 1;
    else
      return N * factorial(N - 1);
  }

  // from the order get the n and m.
  // order = 0.5 * ( n*(n+2) + m )
  // assume the m = 0 and get the n,
  // there will be a unique pair in theory
  std::pair<int, int> get_n_and_l(const std::size_t& order) const {
    int j = static_cast<int>(order);
    int m = 0, n, l, sign_change = 1;
    while (m <= j) {
      l = m * sign_change;
      // get n
      n = static_cast<int>(-1. + std::sqrt(2 * order + 1 - l));
      // check if selected and m and n are correct or not
      if (j == 0.5 * (n * (n + 2) + l)) {
        return {n, l};
        break;
      }

      // incement for the new m
      if (sign_change > 0) {
        m++;
        sign_change = -1;
      } else {
        sign_change = 1;
      }
    }
    // in theory, the n and m should be found, but
    // if didn't get the n and m, then give the fatal error
    if (m == j + 1) {
      fatal_error("the n and m for zernike is not found.");
    }
    return {0, 0};
  }
};

#endif