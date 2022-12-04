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
#ifndef ENTROPY_H
#define ENTROPY_H

#include <array>
#include <cmath>
#include <utils/mpi.hpp>
#include <utils/position.hpp>
#include <vector>

class Entropy {
 public:
  enum Sign { Positive, Negative, Total };

 public:
  Entropy(Position low, Position up, std::array<uint32_t, 3> shp, Sign sgn)
      : lower_corner{low},
        upper_corner{up},
        shape{shp},
        nbins{0},
        dx{},
        dy{},
        dz{},
        bins{},
        total_weight{},
        sign{sgn} {
    nbins = shape[0] * shape[1] * shape[2];
    bins.resize(nbins, 0.);
    total_weight = 0.;
    dx = (upper_corner.x() - lower_corner.x()) / static_cast<double>(shape[0]);
    dy = (upper_corner.y() - lower_corner.y()) / static_cast<double>(shape[1]);
    dz = (upper_corner.z() - lower_corner.z()) / static_cast<double>(shape[2]);
  }
  ~Entropy() = default;

  void add_point(const Position& r, const double& w) {
    // Get bin idecies
    int32_t nx =
        static_cast<int32_t>(std::floor((r.x() - lower_corner.x()) / dx));
    int32_t ny =
        static_cast<int32_t>(std::floor((r.y() - lower_corner.y()) / dy));
    int32_t nz =
        static_cast<int32_t>(std::floor((r.z() - lower_corner.z()) / dz));

    // Make sure position is in a valid bin
    if (nx >= 0 && nx < static_cast<int32_t>(shape[0]) && ny >= 0 &&
        ny < static_cast<int32_t>(shape[1]) && nz >= 0 &&
        nz < static_cast<int32_t>(shape[2])) {
      // Add to total weight, no matter sign being tallied
      this->total_weight += w;

      // Get sign of particle
      Sign p_sign;
      if (w > 0.)
        p_sign = Sign::Positive;
      else
        p_sign = Sign::Negative;

      if (sign == Sign::Total) p_sign = Sign::Total;

      // Add weight to bin if sign agrees
      if (p_sign == sign) {
        this->bins[static_cast<std::size_t>(
            (shape[1] * shape[2]) * static_cast<uint32_t>(nx) +
            (shape[2]) * static_cast<uint32_t>(ny) +
            static_cast<uint32_t>(nz))] += w;
      }

    } else
      std::cout << " Missing entropy particle at " << r << "\n";
  }

  void synchronize_entropy_across_nodes() {
#ifdef ABEILLE_USE_MPI
    // All worker threads must send their generation score to the master.
    // Master must recieve all generations scores from workers and add
    // them to it's own generation score. We do this with MPI_Reduce.
    // We must sadly allocate a recieving buffer to do this however.
    std::vector<double> receiving_buffer;
    double* receiving_buffer_ptr = nullptr;

    // If we are master, we allocate the space in the receiving buffer and
    // set the receiving pointer.
    if (mpi::rank == 0) {
      receiving_buffer.resize(bins.size(), 0.);
      receiving_buffer_ptr = &receiving_buffer[0];
    }

    // Perform reduction
    mpi::synchronize();
    int err = 0;
    err = MPI_Reduce(&bins[0],  // All processes send the tally_gen data
                     receiving_buffer_ptr,  // If we are master, we save sum to
                                            // receiving_buffer
                     static_cast<int>(bins.size()),  // Number of elements
                     mpi::Double,                    // We are sending doubles
                     mpi::Sum,                       // Sum reults
                     0,                              // Results go to master
                     mpi::com);
    mpi::check_error(err, __FILE__, __LINE__);

    // Place results in the receiving_buffer into the tally_gen
    if (mpi::rank == 0) {
      for (std::size_t i = 0; i < bins.size(); i++) {
        bins[i] = receiving_buffer[i];
      }
    }

    // Clear receiving buffer
    receiving_buffer.clear();
    receiving_buffer.shrink_to_fit();

    // Make sure we all have the correct total weight
    double total_weight_tmp = total_weight;
    err = MPI_Reduce(&total_weight_tmp, &total_weight, 1, mpi::Double, mpi::Sum,
                     0, mpi::com);
    mpi::check_error(err, __FILE__, __LINE__);
#endif
  }

  double calculate_entropy() const {
    // Set sum to zero
    double sum = 0.;

    /*for(const auto& b : bins) {
      double p = std::fabs(b) / total_weight;
      if(p != 0.)
        sum -= p*std::log2(p);
    }*/

    for (size_t i = 0; i < this->bins.size(); i++) {
      double p = std::fabs(this->bins[i]) / this->total_weight;

      if (p > 1.0) {
        std::cout << " Negative entropy: p = " << p
                  << ", bin = " << this->bins[i]
                  << ", total_weight = " << this->total_weight << "\n";
      }

      if (p != 0.) {
        sum -= p * std::log2(p);
      }
    }

    return sum;
  }

  void zero() {
    // Zero all bins
    for (auto& b : this->bins) b = 0.;

    this->total_weight = 0.;
  }

  double total_wgt() const { return total_weight; }

 private:
  Position lower_corner;
  Position upper_corner;
  std::array<uint32_t, 3> shape;
  uint32_t nbins;
  double dx, dy, dz;
  std::vector<double> bins;
  double total_weight;
  Sign sign;

};  // Etnropy

#endif  // MG_ENTROPY_H
