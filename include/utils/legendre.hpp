#ifndef LEGENDRE_H
#define LEGENDRE_H

// The following legendre function is required as the standard library in the Mac-OS don't have the Legendre Function. 

#include <cstdint>
#include <utility>

inline double lengendre_orthonormalization(const std::size_t& order,
                                           const double& xmin = -1.0,
                                           const double& xmax = 1.0) {
  return (2.0 * static_cast<double>(order) + 1) / (xmax - xmin);
}

inline double legendre(const std::size_t& order, const double& x) {
  switch (order) {
    case 0:
      return 1.0;
      break;
    case 1:
      return x;
      break;
    case 2:
      return (0.5 * (3.0 * x * x - 1.0));
      break;
    case 3:
      return (0.5 * (5.0 * x * x * x - 3.0 * x));
      break;
    case 4:
      return (0.125 * (35.0 * x * x * x * x - 30.0 * x * x + 3.0));
      break;
    case 5: {
      const double x3 = x * x * x;
      const double x5 = x3 * x * x;
      return (0.125 * (63.0 * x5 - 70.0 * x3 + 15.0 * x));
      break;
    }
    case 6: {
      const double x4 = x * x * x * x;
      const double x6 = x4 * x * x;
      return (0.0625 * (231.0 * x6 - 315.0 * x4 + 105.0 * x * x - 5.0));
      break;
    }
    case 7: {
      const double x3 = x * x * x;
      const double x5 = x3 * x * x;
      const double x7 = x5 * x * x;
      return (0.0625 * (429.0 * x7 - 693.0 * x5 + 315.0 * x3 - 35.0 * x));
      break;
    }
    default: {
      // For 8th or more than 8th order legendre polynomial
      const double x4 = x * x * x * x;
      const double x6 = x4 * x * x;
      double p6 = (0.0625 * (231.0 * x6 - 315.0 * x4 + 105.0 * x * x - 5.0));

      const double x3 = x * x * x;
      const double x5 = x3 * x * x;
      const double x7 = x5 * x * x;
      double p7 = (0.0625 * (429.0 * x7 - 693.0 * x5 + 315.0 * x3 - 35.0 * x));

      unsigned n_ = 7;
      while (n_ < order) {
        p6 = ((2.0 * n_ + 1.0) * x * p7 - n_ * p6) / (n_ + 1.0);
        std::swap(p6, p7);
        n_++;
      }
      return p7;
      break;
    }
  }
}

#endif