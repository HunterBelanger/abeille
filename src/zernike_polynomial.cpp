#include <tallies/zernike_polynomial.hpp>
#include <utils/constants.hpp>

ZernikePolynomials::ZernikePolynomials(std::size_t order)
    : order_(order), max_n_(), Zr_coefficients_(), m_(), n_(), Zr_types_() {
  // Evaluate the polynomial coefficients up to given order
  // Zr_coefficients_ will have the coeffs. stating from lowest power to
  // higher for each order, it will be useful this way as evaluating the power
  // can easily be done.

  // get a maximum size to store the coefficients
  std::pair<std::size_t, int> high_n_and_l = get_n_and_l(order_);
  // store the n for highest order
  // this is required as we have to evaluate upto this power
  max_n_ = high_n_and_l.first;
  std::size_t max_possible_size = 1;
  if (order_ != 0) {
    if (max_n_ % 2 == 0) {
      std::size_t i = (max_n_) / 2;
      max_possible_size = i * (i + 1) * (2 * i + 1) / 3 + i * (i + 1) +
                          (max_n_ + 2) * (max_n_ + 4) / 4;
    } else {
      std::size_t i = (max_n_ + 1) / 2;
      max_possible_size = i * (i + 1) * (2 * i + 1) / 3 + i * (i + 1);
    }
  }
  Zr_coefficients_.reserve(max_possible_size);
  m_.reserve(order_ + 1);
  n_.reserve(order_ + 1);
  Zr_types_.reserve(order_ + 1);

  int l;
  std::size_t n, m;
  for (std::size_t it = 0; it <= order_; it++) {
    std::pair<std::size_t, int> n_and_l = get_n_and_l(it);
    n = n_and_l.first;
    l = n_and_l.second;
    m = static_cast<std::size_t>(std::abs(l));
    ZernikeType type;
    if (l < 0) {
      type = ZernikeType::Odd;
    } else {
      type = ZernikeType::Even;
    }
    m_.push_back(m);
    n_.push_back(n);
    Zr_types_.push_back(type);

    std::size_t max_k = static_cast<std::size_t>((n - m) / 2);
    std::vector<double> single_order_Zr_coeff;
    single_order_Zr_coeff.reserve(max_k + 1);
    double sign_change = 1.;
    if (static_cast<int>(max_k % 2) == 1) sign_change = -1.;

    std::size_t k;
    for (std::size_t i = 0; i <= max_k; i++) {
      k = max_k - i;
      const double numeriator = factorial(n - k);
      const std::size_t avg_n_m = (n + m) / 2;
      const std::size_t mid_n_m = (n - m) / 2;
      const double denominator =
          factorial(k) * factorial(avg_n_m - k) * factorial(mid_n_m - k);

      const double coeff = sign_change * numeriator / denominator;
      Zr_coefficients_.push_back(coeff);
      sign_change = -1 * sign_change;
    }
  }
  Zr_coefficients_.shrink_to_fit();
}


std::vector<double> ZernikePolynomials::evaluate_zernikes(
    const double r, const double theta) const {
  // evaluate the all the powers up maximum order, i.e, max_n_
  std::vector<double> r_powers;
  r_powers.reserve(max_n_ + 1);
  double r_k = 1.0;
  r_powers.push_back(r_k);
  for (std::size_t i = 1; i <= max_n_; i++) {
    r_k *= r;
    r_powers.push_back(r_k);
  }

  // evaluate the each order polynomial value and store in the vector
  std::size_t it_coeff = 0;
  std::vector<double> zr_values;
  zr_values.reserve(order_ + 1);
  for (std::size_t order = 0; order <= order_; order++) {
    const std::size_t m = m_[order];
    const std::size_t n = n_[order];
    const std::size_t order_coeff_size = (n - m) / 2 + 1;
    const ZernikeType zr_type = Zr_types_[order];
    double value = 0.0;
    std::size_t k = m;
    for (std::size_t i = 0; i < order_coeff_size; i++) {
      value += Zr_coefficients_[it_coeff] * r_powers[k];
      k += 2;
      it_coeff++;
    }
    const double mtheta = static_cast<double>(m) * theta;
    if (zr_type == ZernikeType::Even) {
      value *= std::cos(mtheta);
    } else {
      value *= std::sin(mtheta);
    }
    zr_values.push_back(value);
  }

  return zr_values;
}

// orthonormalsation constnat will be achieved by the inverse of square of L-2
// norm.
double ZernikePolynomials::orthonormalization_constant(
    const std::size_t& order) const {
  std::pair<int, int> n_and_l = get_n_and_l(order);
  const double n = static_cast<double>(n_and_l.first);
  const double l = static_cast<double>(n_and_l.second);
  if (l == 0) {
    return 1. / (n + 1.);
  }
  return 0.5 / (n + 1.);
}