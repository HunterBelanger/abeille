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

  std::vector<double> evaluate_zernikes(const double r,
                                        const double theta) const;

  double orthonormalization_constant(const std::size_t& order) const;

 private:
  std::size_t order_, max_n_;
  std::vector<double> Zr_coefficients_;
  std::vector<std::size_t> m_;
  std::vector<std::size_t> n_;
  std::vector<ZernikeType> Zr_types_;

  // function for the factorial
  double factorial(std::size_t N) const {
    if (N == 1 || N == 0)
      return 1;
    else
      return static_cast<double>(N) * factorial(N - 1);
  }

  // from the order get the n and l.
  // order = 0.5 * ( n*(n+2) + l )
  // assume the n = 0 and get the l, if not satisfied then increase the l
  // there will be a unique pair in theory
  std::pair<std::size_t, int> get_n_and_l(const std::size_t& order) const {
    int n = 0;
    while (n <= static_cast<int>(order)) {
      int l = 2 * static_cast<int>(order) - n * (n + 2);
      if (std::abs(l) <= n) return {static_cast<std::size_t>(n), l};
      n++;
    }
    return {0, 0};
  }
};

#endif