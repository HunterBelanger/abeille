/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *============================================================================*/
#include <simulation/collision_mesh_tally.hpp>
#include <simulation/mesh_tally.hpp>
#include <simulation/track_length_mesh_tally.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <vector>

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
      net_weight(1.),
      energy_bounds(ebounds),
      fname(fname),
      tally_gen(),
      tally_avg(),
      tally_var() {
  dx = (r_hi.x() - r_low.x()) / static_cast<double>(Nx);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(Ny);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(Nz);

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

  // All worker threads must send their generation score to the master.
  // Master must recieve all generations scores from workers and add
  // them to it's own generation score.

  mpi::Reduce_sum(&tally_gen[0], static_cast<int>(tally_gen.size()), 0);

  // Only try to update average and variance is we are master, as worker
  // processes don't have copies of this data, so it will seg-fault.
  if (mpi::rank == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < tally_gen.size(); i++) {
      // Get new average
      double old_avg = tally_avg[i];
      double val = tally_gen[i] * multiplier;
      double avg = old_avg + (val - old_avg) / dg;
      tally_avg[i] = avg;

      // Get new variance
      double var = tally_var[i];
      var = var + (((val - old_avg) * (val - avg) - (var)) / dg);
      tally_var[i] = var;
    }
  }
}

void MeshTally::clear_generation() { tally_gen.fill(0.); }

void MeshTally::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank == 0) {
    std::string mssg = " Writing " + fname + " tally file...\n";
    Output::instance()->write(mssg);

    std::ofstream file(fname + "_mesh.txt");

    // First write coordinates and number of groups
    file << " X:";
    for (uint64_t i = 0; i <= Nx; i++) {
      double x = (static_cast<double>(i) * dx) + r_low.x();
      file << x;
      if (i != Nx)
        file << ",";
      else
        file << "\n";
    }

    file << " Y:";
    for (uint64_t i = 0; i <= Ny; i++) {
      double y = (static_cast<double>(i) * dy) + r_low.y();
      file << y;
      if (i != Ny)
        file << ",";
      else
        file << "\n";
    }

    file << " Z:";
    for (uint64_t i = 0; i <= Nz; i++) {
      double z = (static_cast<double>(i) * dz) + r_low.z();
      file << z;
      if (i != Nz)
        file << ",";
      else
        file << "\n";
    }

    file << " ENERGY:";
    for (std::size_t i = 0; i < energy_bounds.size(); i++) {
      file << energy_bounds[i];
      if (i != energy_bounds.size() - 1)
        file << ",";
      else
        file << "\n";
    }

    file << " QUANTITY: " << this->quantity_str() << "\n";

    file.close();

    // Convert flux_var to the error on the mean
    for (size_t l = 0; l < tally_var.size(); l++)
      tally_var[l] = std::sqrt(tally_var[l] / static_cast<double>(g));

    tally_avg.save(fname + "_avg.npy");
    tally_var.save(fname + "_err.npy");
  }
}
