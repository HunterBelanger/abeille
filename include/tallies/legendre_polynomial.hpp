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
  std::vector<double> legendre_coeff_;
  std::size_t order_;

  // function to evaluate the factorial calculations for coefficients for given order and term
  double coeff_factorial_evaluation(const std::size_t& n, const std::size_t& k){
    double value = 1.; 
    std::size_t Nk = n-k;
    for (std::size_t i = 1; i <= Nk; i++ ){
      double factorial_n_2k = 1.;
      double factorial_k = 1.;
      if ( i < k + 1 ){
        factorial_k = i;
      }
      if ( (n-2*k) + 1 > i ){
        factorial_n_2k = i;
      }
      value *= (n-k  + i) / ( factorial_k * factorial_n_2k );
    }
    return value;
  }
};

#endif