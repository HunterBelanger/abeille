#ifndef LEGENDRE_POLYNOMIALS_H
#define LEGENDRE_POLYNOMIALS_H

#include <vector>

//===================================
// Legendre Polynomials from order 0 to N order.
//
class LegendrePolynomials {
 public:
 
  LegendrePolynomials(std::size_t order);

  std::vector<double> evaluate_legendres(const double x) const;

 private:
  std::vector<double> legendre_coeff_;
  std::size_t order_;

  // function to evaluate the factorial calculations for coefficients for given
  // order and term
  double coeff_factorial_evaluation(const std::size_t& n,
                                    const std::size_t& k) const {
    double value = 1.;
    std::size_t Nk = n - k;
    for (std::size_t i = 1; i <= Nk; i++) {
      if (i < k + 1) {
        value /= static_cast<double>(i);
      }
      if (i < (n - 2 * k) + 1) {
        value /= static_cast<double>(i);
      }
      value *= static_cast<double>(n - k + i);
    }
    return value;
  }
};

#endif