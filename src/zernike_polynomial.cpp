#include <tallies/zernike_polynomial.hpp>
#include <utils/constants.hpp>

ZernikePolynomials::ZernikePolynomials(std::size_t order)
    : order_(order), max_n_(), Zr_coefficients_(), m_(), n_(), Zr_types_() {
  // Evaluate the polynomial coefficients up to given order
  // Zr_coefficients_ will be having the coeff. stating from lowest power to
  // higher for each order, it will be useful this way as evaluating the power
  // can easily be done.

  // get a maximum size to store the coefficients
  std::pair<int, int> high_n_and_l = get_n_and_l(order_);
  // store the n for highest order
  // this is required as we have to evaluate upto this power
  max_n_ = static_cast<std::size_t>(high_n_and_l.first);
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

  int m, n, l;
  for (std::size_t it = 0; it <= order_; it++) {
    std::pair<int, int> n_and_l = get_n_and_l(it);
    n = n_and_l.first;
    l = n_and_l.second;
    m = std::abs(l);
    ZernikeType type;
    if (l < 0) {
      type = ZernikeType::Odd;
    } else {
      type = ZernikeType::Even;
    }
    m_.push_back(static_cast<std::size_t>(m));
    n_.push_back(static_cast<std::size_t>(n));
    Zr_types_.push_back(type);

    std::size_t max_k = static_cast<std::size_t>((n - m) * 0.5);
    std::vector<double> single_order_Zr_coeff;
    single_order_Zr_coeff.reserve(max_k + 1);
    double sign_change = 1.;

    if (static_cast<int>(max_k % 2) == 1) sign_change = -1.;

    std::size_t k;
    for (std::size_t i = 0; i <= max_k; i++) {
      k = max_k - i;
      const double numeriator = factorial(n - k);
      const std::size_t avg_n_m = static_cast<std::size_t>((n + m) * 0.5);
      const std::size_t mid_n_m = static_cast<std::size_t>((n - m) * 0.5);
      const double denominator =
          factorial(k) * factorial(avg_n_m - k) * factorial(mid_n_m - k);

      const double coeff = sign_change * numeriator / denominator;
      Zr_coefficients_.push_back(coeff);
      sign_change = -1 * sign_change;
    }
  }
  Zr_coefficients_.shrink_to_fit();
}

double ZernikePolynomials::evaluate_zernike_at_order(
    const double x, const double theta, const std::size_t order) const {
  if (order_ < order) {
    fatal_error(
        "the providede order for the zernike-polynomial is out of the stored "
        "max-order.");
  }

  // get the start and end location of desired order's coefficients
  std::size_t it_coeff_start = 0;
  std::size_t it_coeff_end = 0;
  if (order != 0) {
    for (std::size_t ord = 0; ord < order; ord++) {
      std::pair<int, int> n_and_l = get_n_and_l(ord);
      const std::size_t n = static_cast<std::size_t>(n_and_l.first);
      const std::size_t m = static_cast<std::size_t>(std::abs(n_and_l.second));
      it_coeff_start += static_cast<std::size_t>((n - m) * 0.5 + 1);
    }
    std::pair<int, int> n_and_l = get_n_and_l(order);
    const std::size_t n = static_cast<std::size_t>(n_and_l.first);
    const std::size_t m = static_cast<std::size_t>(std::abs(n_and_l.second));
    it_coeff_end = it_coeff_start + static_cast<std::size_t>((n - m) * 0.5 + 1);
  }

  const std::size_t m = m_[order];
  double x_k = std::pow(x, m);
  const double x_k2 = x * x;

  double value = 0.;
  for (std::size_t i = it_coeff_start; i < it_coeff_end; i++) {
    value += Zr_coefficients_[i] * x_k;
    x_k *= x_k2;
  }
  double mtheta = static_cast<double>(m) * theta;
  if (Zr_types_[order] == ZernikeType::Even) {
    value *= std::cos(mtheta);
  } else {
    value *= std::sin(mtheta);
  }

  return value;
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
    const std::size_t order_coeff_size =
        static_cast<std::size_t>((n - m) * 0.5 + 1.);
    const ZernikeType zr_type = Zr_types_[order];
    double value = 0.0;
    std::size_t k = m;
    for (std::size_t i = 0; i < order_coeff_size; i++) {
      value += Zr_coefficients_[it_coeff] * r_powers[k];
      k += 2;
      it_coeff++;
    }
    double mtheta = static_cast<double>(m) * theta;
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