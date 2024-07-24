#ifndef LEGENDRE_POLYNOMIALS_H
#define LEGENDRE_POLYNOMIALS_H

#include <vector>

//===================================
// Legendre Polynomials from order 0 to N order.
//
class LegendrePolynomials {
 public:
  LegendrePolynomials() = default;
  LegendrePolynomials(std::size_t order);

  //~LegendrePolynomials() = default;

  std::vector<double> evaluate_legendres(const double x) const;

 private:
  std::vector<std::vector<double>> legendre_coeff_;
  std::size_t order_;

  // function for the factorial
  std::size_t factorial(std::size_t N) {
    if (N == 1 || N == 0)
      return 1;
    else
      return N * factorial(N - 1);
  }
};

#endif