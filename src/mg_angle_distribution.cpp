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
#include <materials/mg_angle_distribution.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

#include <algorithm>
#include <sstream>

#include <PapillonNDL/linearize.hpp>
#include <PapillonNDL/pctable.hpp>
#include <PapillonNDL/tabulated_1d.hpp>

MGAngleDistribution::MGAngleDistribution()
    : mu_({-1., 1.}),
      pdf_({0.5, 0.5}),
      cdf_({0., 1.}),
      abs_pdf_(mu_, pdf_, cdf_, pndl::Interpolation::LinLin) {}

MGAngleDistribution::MGAngleDistribution(const std::vector<double>& mu,
                                         const std::vector<double>& pdf,
                                         const std::vector<double>& cdf)
    : mu_(mu),
      pdf_(pdf),
      cdf_(cdf),
      abs_pdf_({-1.0, 1.0}, {0.5, 0.5}, {0.0, 1.0},
               pndl::Interpolation::LinLin) {
  // Make sure good mu bounds
  if (mu_.front() < -1.) {
    fatal_error("Angle limit less than -1.");
  }

  if (mu_.back() > 1.) {
    fatal_error("Angle limit greater than 1.");
  }

  // Make sure mu is sorted
  if (std::is_sorted(mu_.begin(), mu_.end()) == false) {
    fatal_error("Mu values are not sorted.");
  }

  // Make sure PDF is positive
  // If portions of the pdf are found to be negative, then a weight modifier
  // will be used to sample the distribuiton with importance sampling.
  for (const auto& p : pdf_) {
    if (p < 0.) {
      warning("Encountered angular distribution with a negative PDF.");
      pdf_is_neg = true;
    }
  }

  // Setup the variables for negative pdf distribution
  if (pdf_is_neg == true) {
    // abs_neg_pdf will store absolute value the negative distribuion
    std::vector<double> abs_neg_pdf_;
    abs_neg_pdf_.reserve(pdf_.size());

    for (const auto& p : pdf_) {
      abs_neg_pdf_.push_back(std::abs(p));
    }

    pndl::Tabulated1D pdf_original(pndl::Interpolation::LinLin, mu_, pdf_);

    // Lambda function to get the values from pdf_orginal, containing the
    // absolute pdf values
    auto abs_pdf_function = [&pdf_original](double x) {
      return std::abs(pdf_original(x));
    };

    pndl::Tabulated1D abs_pdf_tabulated_ =
        pndl::linearize(mu_, abs_neg_pdf_, abs_pdf_function);

    // area under the abs distribution
    abs_weight_mod_ = abs_pdf_tabulated_.integrate(mu_.front(), mu_.back());
    const double inverse_abs_pdf_area = 1. / abs_weight_mod_;

    // pdf reformation based on the abs_pdf 
    abs_neg_pdf_ = abs_pdf_tabulated_.y();
    abs_neg_pdf_[0] *= inverse_abs_pdf_area;
    
    // cdf corresponds to abs negative distribution
    std::vector<double> abs_neg_cdf_(abs_neg_pdf_.size(), 0.);
    
    for (std::size_t i = 1; i < abs_pdf_tabulated_.x().size(); i++) {
      abs_neg_pdf_[i] *= inverse_abs_pdf_area;  // pdf normalization
      abs_neg_cdf_[i] =
          abs_neg_cdf_[i - 1] +
          0.5 * (abs_neg_pdf_[i] + abs_neg_pdf_[i - 1]) *
              (abs_pdf_tabulated_.x()[i] - abs_pdf_tabulated_.x()[i - 1]);
    }

    abs_pdf_ = pndl::PCTable(abs_pdf_tabulated_.x(), abs_neg_pdf_, abs_neg_cdf_,
                             pndl::Interpolation::LinLin);
  } else
    // Make sure CDF is sorted and > 0
    if (std::is_sorted(cdf_.begin(), cdf_.end()) == false) {
      fatal_error("CDF is not sorted.");
    }
      
    if (cdf_.front() < 0.) {
      fatal_error("CDF is less than 0.");
    }
  }

  // Make sure CDF starts at 0, and ends at 1
  if (cdf_.front() != 0.) {
    fatal_error("First CDF value is not 0.");
  }

  if (cdf_.back() != 1.) {
    fatal_error("Last CDF value is not 1.");
  }
}
