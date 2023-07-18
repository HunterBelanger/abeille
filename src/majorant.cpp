/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
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
#include <materials/ce_nuclide.hpp>
#include <materials/material_helper.hpp>
#include <materials/nuclide.hpp>
#include <utils/majorant.hpp>
#include <utils/settings.hpp>

#include <boost/unordered/unordered_flat_map.hpp>

#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/st_coherent_elastic.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>
#include <PapillonNDL/st_thermal_scattering_law.hpp>

#include <algorithm>
#include <cmath>

std::pair<std::vector<double>, std::vector<double>> make_majorant_xs() {
  // Must first create a unionized energy grid. How this is done depends on
  // whether or not we are in continuous energy or multi-group mode.
  if (settings::energy_mode == settings::EnergyMode::CE) {
    // First, we must construct a unionized energy grid, for all nuclides in the
    // problem. We do this by iterating through all nuclides. In CE mode, we
    // should be able to safely case a Nuclide pointer to a CENuclide pointer.
    std::vector<double> union_energy_grid;

    // Get a map of all URR random values set to 1, for finding the majorant
    boost::unordered_flat_map<uint32_t, std::optional<double>> urr_rands;
    for (const auto& za : zaids_with_urr) urr_rands[za] = 1.;

    // Iterate through all nuclides
    for (const auto& id_nucld_pair : nuclides) {
      Nuclide* nuc = id_nucld_pair.second.get();

      // If we are in CE mode, this static cast should be safe.
      CENuclide* cenuc = static_cast<CENuclide*>(nuc);

      // Get the nuclide's energy grid
      std::vector<double> cenuc_grid = cenuc->cedata()->energy_grid().grid();
      if (cenuc->tsl()) {
        // If we have a TSL, we need to change cenuc_grid to add those points.
        // First, get ref to IncoherentInelastic
        const pndl::STIncoherentInelastic& II =
            cenuc->tsl()->incoherent_inelastic();
        const pndl::STCoherentElastic& CE = cenuc->tsl()->coherent_elastic();
        std::size_t NEtsl = II.xs().x().size() + (2 * CE.bragg_edges().size());
        cenuc_grid.reserve(cenuc_grid.size() + NEtsl);

        // First, add all II points
        for (std::size_t i = 0; i < II.xs().x().size(); i++) {
          cenuc_grid.push_back(II.xs().x()[i]);
        }

        // Now add all CE bragg edges
        for (std::size_t i = 0; i < CE.bragg_edges().size(); i++) {
          cenuc_grid.push_back(std::nextafter(CE.bragg_edges()[i], 0.));
          cenuc_grid.push_back(CE.bragg_edges()[i]);
        }

        // Now sort the grid
        std::sort(cenuc_grid.begin(), cenuc_grid.end());
      }

      // Get the URR grid points
      if (cenuc->has_urr()) {
        cenuc_grid.reserve(cenuc_grid.size() + cenuc->urr_energy_grid().size());

        // Get the URR energy points
        for (const auto& urr_E : cenuc->urr_energy_grid()) {
          cenuc_grid.push_back(urr_E);
        }

        // Now sort the grid
        std::sort(cenuc_grid.begin(), cenuc_grid.end());
      }

      // Now we perform the union operation
      std::vector<double> tmp;
      std::set_union(union_energy_grid.begin(), union_energy_grid.end(),
                     cenuc_grid.begin(), cenuc_grid.end(),
                     std::back_inserter(tmp));
      union_energy_grid = tmp;

      // We now throw the grid into a set temporarily, to remove the duplicates
      std::set<double> tmp_union_energy_grid(union_energy_grid.begin(),
                                             union_energy_grid.end());
      union_energy_grid.assign(tmp_union_energy_grid.begin(),
                               tmp_union_energy_grid.end());
    }

    // With the energy grid unionized, we can go about finding majorants !
    std::vector<double> maj_xs(union_energy_grid.size(), 0.);
    for (const auto& material : materials) {
      // Make a MaterialHelper, to evaluate the material XS for us
      MaterialHelper mat(material.second.get(), union_energy_grid.front());
      mat.set_urr_rand_vals(urr_rands);

      // Now iterate through all energy points. If xs is larger, update value
      for (std::size_t i = 0; i < union_energy_grid.size(); i++) {
        const double Ei = union_energy_grid[i];
        const double xsi = mat.Et(Ei);

        if (xsi > maj_xs[i]) maj_xs[i] = xsi;
      }
    }

    // Multiply all majorant values by a small safety factor
    for (auto& xsmaj : maj_xs) xsmaj *= 1.01;

    return {union_energy_grid, maj_xs};
  } else {
    // We are in multi-group mode. Here, the energy-bounds are kept in the
    // settings, so we can construct something with that

    std::vector<double> egrid;
    egrid.push_back(settings::energy_bounds[0]);
    if (egrid.front() == 0.) egrid.front() = 1.E-11;

    for (size_t i = 1; i < settings::energy_bounds.size() - 1; i++) {
      egrid.push_back(settings::energy_bounds[i]);
      egrid.push_back(settings::energy_bounds[i]);
    }

    egrid.push_back(settings::energy_bounds.back());

    // This now has created the vector egrid which will look something like this
    // for the case of 5 energy groups.
    // [0., 1.,   1., 2.,   2., 3.,   3., 4.,   4., 5.]
    // This works, because the energy of multi-group particles should always be
    // inbetween the bounds for the group.

    // Now we need to make a vector which will contian the majorant cross
    // cross sections for each group.
    std::vector<double> maj_xs(egrid.size(), 0.);

    // We loop through materials
    for (const auto& material : materials) {
      // Then we loop through energies
      MaterialHelper mat(material.second.get(), 1.);

      for (uint32_t g = 0; g < settings::ngroups; g++) {
        // Get the energy at the mid-point for the group
        size_t i = g * 2;
        double Eg = 0.5 * (egrid[i] + egrid[i + 1]);

        double xs = mat.Et(Eg);
        if (xs > maj_xs[i]) {
          maj_xs[i] = xs;
          maj_xs[i + 1] = xs;
        }
      }
    }

    return {egrid, maj_xs};
  }
}