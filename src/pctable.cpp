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
#include <utils/error.hpp>
#include <utils/pctable.hpp>

#include <vector>

pndl::PCTable make_pctable(const YAML::Node& node) {
  // Get Values
  if (!node["values"] || !node["values"].IsSequence()) {
    fatal_error("No values list in tabulated distribution.");
  }
  std::vector<double> values = node["values"].as<std::vector<double>>();

  if (values.size() < 2) {
    fatal_error(
        "Must have at least 2 tabulated points for tabulated distribution.");
  }

  if (!std::is_sorted(values.begin(), values.end())) {
    fatal_error("All values in tabulated distrubtion must be sorted.");
  }

  // Get PDF
  if (!node["pdf"] || !node["pdf"].IsSequence()) {
    fatal_error("No pdf list in tabulated distribution.");
  }
  std::vector<double> pdf = node["pdf"].as<std::vector<double>>();

  if (values.size() != pdf.size()) {
    fatal_error("Values and pdf lists must have the same size.");
  }

  for (const auto& p : pdf) {
    if (p < 0.) {
      fatal_error("All pdf entries in tabulated distrubtion must be >= 0.");
    }
  }

  // Construct CDF
  std::vector<double> cdf(values.size(), 0.);
  for (std::size_t i = 0; i < cdf.size() - 1; i++) {
    cdf[i + 1] =
        cdf[i] + 0.5 * (values[i + 1] - values[i]) * (pdf[i + 1] + pdf[i]);
  }

  // Normalize PDF and CDF
  const double cdf_max = cdf.back();
  for (std::size_t i = 0; i < cdf.size(); i++) {
    pdf[i] /= cdf_max;
    cdf[i] /= cdf_max;
  }

  return pndl::PCTable(values, pdf, cdf, pndl::Interpolation::LinLin);
}