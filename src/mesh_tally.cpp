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
#include <tallies/collision_mesh_tally.hpp>
#include <tallies/mesh_tally.hpp>
#include <tallies/track_length_mesh_tally.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <tallies/general_tally.hpp>

#include <set>
#include <vector>

const static std::set<std::string> disallowed_tally_names{
    "families",
    "pair-dist-sqrd",
    "entropy",
    "total-pre-cancel-entropy",
    "neg-pre-cancel-entropy",
    "pos-pre-cancel-entropy",
    "total-post-cancel-entropy",
    "neg-post-cancel-entropy",
    "pos-post-cancel-entropy",
    "empty-entropy-frac",
    "Nnet",
    "Ntot",
    "Npos",
    "Nneg",
    "Wnet",
    "Wtot",
    "Wpos",
    "Wneg",
    "kcol",
    "ktrk",
    "kabs",
    "leakage",
    "mig-area"};

MeshTally::MeshTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                     uint64_t nz, const std::vector<double>& ebounds,
                     std::string fname)
    : r_low{low},
      r_hi{hi},
      Nx{nx},
      Ny{ny},
      Nz{nz},
      g(0),
      dx(),
      dy(),
      dz(),
      dx_inv(),
      dy_inv(),
      dz_inv(),
      net_weight(1.),
      energy_bounds(ebounds),
      fname(fname),
      tally_gen(),
      tally_avg(),
      tally_var() {
  // Make sure the name is allowed.
  if (disallowed_tally_names.contains(this->fname)) {
    fatal_error("The tally name " + this->fname + " is reserved.");
  }

  dx = (r_hi.x() - r_low.x()) / static_cast<double>(Nx);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(Ny);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(Nz);
  dx_inv = 1. / dx;
  dy_inv = 1. / dy;
  dz_inv = 1. / dz;

  uint32_t Ne = static_cast<uint32_t>(energy_bounds.size() - 1);

  // Allocate and fill arrays to zero
  tally_gen.reallocate({Ne, Nx, Ny, Nz});
  tally_gen.fill(0.);

  // Only allocate average and variance if we are the master !
  if (mpi::rank == 0) {
    tally_avg.reallocate({Ne, Nx, Ny, Nz});
    tally_avg.fill(0.);

    tally_var.reallocate({Ne, Nx, Ny, Nz});
    tally_var.fill(0.);
  }
}

void MeshTally::set_net_weight(double W) { net_weight = W; }

void MeshTally::record_generation(double multiplier) {
  // Advance the number of generations
  g++;

  const double dg = static_cast<double>(g);
  const double invs_dg = 1. / dg;

  // All worker threads must send their generation score to the master.
  // Master must recieve all generations scores from workers and add
  // them to it's own generation score.
  mpi::Reduce_sum(tally_gen.data_vector(), 0);

  // Only try to update average and variance is we are master, as worker
  // processes don't have copies of this data, so it will seg-fault.
  if (mpi::rank == 0) {
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < tally_gen.size(); i++) {
      // Get new average
      double old_avg = tally_avg[i];
      double val = tally_gen[i] * multiplier;
      double avg = old_avg + (val - old_avg) * invs_dg;
      tally_avg[i] = avg;

      // Get new variance
      double var = tally_var[i];
      var = var + ((val - old_avg) * (val - avg) - (var)) * invs_dg;
      tally_var[i] = var;

      // std::cout<<"Given :"<<tally_avg[i]<<"\n";
    }
  }
}

void MeshTally::clear_generation() { tally_gen.fill(0.); }

void MeshTally::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank != 0) return;

  auto& h5 = Output::instance().h5();

  // Create the group for the tally
  auto tally_grp = h5.createGroup("results/" + this->fname);

  // First write coordinates and number of groups
  std::vector<double> x_bounds(Nx + 1, 0.);
  for (std::size_t i = 0; i <= Nx; i++) {
    x_bounds[i] = (static_cast<double>(i) * dx) + r_low.x();
  }
  tally_grp.createAttribute("x-bounds", x_bounds);

  std::vector<double> y_bounds(Ny + 1, 0.);
  for (std::size_t i = 0; i <= Ny; i++) {
    y_bounds[i] = (static_cast<double>(i) * dy) + r_low.y();
  }
  tally_grp.createAttribute("y-bounds", y_bounds);

  std::vector<double> z_bounds(Nz + 1, 0.);
  for (std::size_t i = 0; i <= Nz; i++) {
    z_bounds[i] = (static_cast<double>(i) * dz) + r_low.z();
  }
  tally_grp.createAttribute("z-bounds", z_bounds);

  tally_grp.createAttribute("energy-bounds", energy_bounds);

  // Save the quantity
  tally_grp.createAttribute("quantity", this->quantity_str());

  if (this->quantity_str() == "mt") {
    tally_grp.createAttribute("mt", this->mt());
  }

  // Save the estimator
  tally_grp.createAttribute("estimator", this->estimator_str());

  // Convert flux_var to the error on the mean
  for (size_t l = 0; l < tally_var.size(); l++)
    tally_var[l] = std::sqrt(tally_var[l] / static_cast<double>(g));

  // Add data sets for the average and the standard deviation
  auto avg_dset =
      tally_grp.createDataSet<double>("avg", H5::DataSpace(tally_avg.shape()));
  avg_dset.write_raw(&tally_avg[0]);

  auto std_dset =
      tally_grp.createDataSet<double>("std", H5::DataSpace(tally_var.shape()));
  std_dset.write_raw(&tally_var[0]);
}
