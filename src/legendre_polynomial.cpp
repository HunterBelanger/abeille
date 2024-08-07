#include <tallies/legendre_polynomial.hpp>

#include <iostream>

LegendrePolynomials::LegendrePolynomials(std::size_t order)
    : legendre_coeff_(), order_(order) {
  // evaluate the coefficient for each and every order
  // from zero to given order

  // get the total number of coefficient to store
  std::size_t coeff_length = 1;
  if (order_ != 0) {
    if (order_ % 2 == 0) {
      const std::size_t length =
          static_cast<std::size_t>((static_cast<double>(order_) - 1.) * 0.5 + 1.);
      coeff_length =
          length * (length + 1) + static_cast<std::size_t>(static_cast<double>(order_) * 0.5 + 1);
    } else {
      const std::size_t length = static_cast<std::size_t>(static_cast<double>(order_) * 0.5 + 1);
      coeff_length = length * (length + 1);
    }
  }
  legendre_coeff_.reserve(coeff_length);

  double inv_2 = 1.;
  std::size_t low_pow = 1;

  for (std::size_t i = 0; i <= order_; i++) {
    if (low_pow == 1) {
      low_pow = 0;
    } else if (low_pow == 0) {
      low_pow = 1;
    }
    std::size_t high_pow = static_cast<std::size_t>(static_cast<double>(i) * 0.5);
    double sign_change = 1.;
    if (high_pow % 2 == 1) sign_change = -1.;

    for (std::size_t j = 0; j <= high_pow; j++) {
      std::size_t k = high_pow - j;
      double value = inv_2 * sign_change * coeff_factorial_evaluation(i, k);
      sign_change *= -1;
      legendre_coeff_.push_back(value);
    }
    inv_2 *= 0.5;
  }

  std::vector<double> sol = evaluate_legendres(0.1);
  for (std::size_t i = 0; i < sol.size(); i++) {
    std::cout << "order = " << i << "\t\t" << sol[i] << "\n";
  }
}

std::vector<double> LegendrePolynomials::evaluate_legendres(
    const double x) const {
  // get the powers up to that order
  std::vector<double> x_powers;
  x_powers.reserve(order_ + 1);
  double x_pow = 1.;
  for (std::size_t i = 0; i <= order_; i++) {
    x_powers.push_back(x_pow);
    x_pow *= x;
  }
  // evaluate the legendre polynomials of each order
  std::vector<double> legndre_values;
  legndre_values.reserve(order_ + 1);
  std::size_t low_pow = 1;
  // iteration over the coeff.
  std::size_t it_coeff = 0;
  for (std::size_t i = 0; i <= order_; i++) {
    if (low_pow == 0)
      low_pow = 1;
    else
      low_pow = 0;

    double sum = 0.;
    for (std::size_t k = low_pow; k <= i; k += 2) {
      sum += legendre_coeff_[it_coeff] * x_powers[k];
      it_coeff++;
    }
    legndre_values.push_back(sum);
  }
  return legndre_values;
}