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
#ifndef MPI_H
#define MPI_H

#include <cstdint>
#include <vector>

#ifdef ABEILLE_USE_MPI
#include <mpi.h>

#include <simulation/particle.hpp>
#endif

#include <utils/timer.hpp>

namespace mpi {

// Timer for all MPI operations
extern Timer timer;

#ifdef ABEILLE_USE_MPI
using Com = MPI_Comm;
using DType = MPI_Datatype;
using OpType = MPI_Op;

extern const Com com;
extern const DType Bool;
extern const DType Int;
extern const DType Double;
extern const DType UInt64;
extern DType BParticle;

extern const OpType Sum;
extern const OpType And;
extern const OpType Or;

//==============================================================================
// MPI Data Types
template <typename T>
DType dtype();

template <>
inline DType dtype<bool>() {
  return Bool;
}

template <>
inline DType dtype<int>() {
  return Int;
}

template <>
inline DType dtype<double>() {
  return Double;
}

template <>
inline DType dtype<uint64_t>() {
  return UInt64;
}

template <>
inline DType dtype<BankedParticle>() {
  return BParticle;
}
#endif

extern std::vector<uint64_t> node_nparticles;
extern std::vector<uint64_t> node_nparticles_noise;

extern int size;
extern int rank;

void initialize_mpi(int* argc, char*** argv);
void finalize_mpi();
void abort_mpi();
void synchronize();
void check_error(int err, const char* file, int line);

//==============================================================================
// MPI Operations
template <typename T>
void Bcast(T& val, int root, const char* file = __FILE__, int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err = MPI_Bcast(&val, 1, dtype<T>(), root, com);
    check_error(err, file, line);
    timer.stop();
  }
#else
  (void)val;
  (void)root;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Allreduce_sum(T& val, const char* file = __FILE__, int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    T tmp_send = val;
    int err = MPI_Allreduce(&tmp_send, &val, 1, dtype<T>(), Sum, com);
    check_error(err, file, line);
    timer.stop();
  }
#else
  (void)val;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Allreduce_or(T& val, const char* file = __FILE__, int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    T tmp_send = val;
    int err = MPI_Allreduce(&tmp_send, &val, 1, dtype<T>(), Or, com);
    check_error(err, file, line);
    timer.stop();
  }
#else
  (void)val;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Reduce_sum(T& val, int root, const char* file = __FILE__,
                int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    T tmp_send = val;
    int err = MPI_Reduce(&tmp_send, &val, 1, dtype<T>(), Sum, root, com);
    check_error(err, file, line);
    timer.stop();
  }
#else
  (void)val;
  (void)root;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Reduce_sum(std::vector<T>& vals, int root, const char* file = __FILE__,
                int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    std::vector<T> tmp_rcv;

    if (rank == root) {
      // Only allocate the receving array if we are the reciever
      tmp_rcv.resize(vals.size());
    }

    int err = MPI_Reduce(&vals[0], &tmp_rcv[0], static_cast<int>(vals.size()),
                         dtype<T>(), Sum, root, com);
    check_error(err, file, line);

    if (rank == root) {
      for (std::size_t i = 0; i < vals; i++) {
        vals[i] = tmp_rcv[i];
      }
    }
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Reduce_sum(T* vals, int count, int root, const char* file = __FILE__,
                int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    std::vector<T> tmp_rcv;

    if (rank == root) {
      // Only allocate the receving array if we are the reciever
      tmp_rcv.resize(static_cast<std::size_t>(count));
    }

    int err = MPI_Reduce(vals, &tmp_rcv[0], count, dtype<T>(), Sum, root, com);
    check_error(err, file, line);

    if (rank == root) {
      for (std::size_t i = 0; i < tmp_rcv.size(); i++) {
        vals[i] = tmp_rcv[i];
      }
    }
    timer.stop();
  }
#else
  (void)vals;
  (void)count;
  (void)root;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Gatherv(std::vector<T>& vals, int root, const char* file = __FILE__,
             int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    // First, we need to know how many things each rank has
    std::vector<int> sizes(static_cast<std::size_t>(size), 0);
    sizes[static_cast<std::size_t>(rank)] = static_cast<int>(vals.size());
    for (std::size_t n = 0; n < static_cast<std::size_t>(size); n++) {
      Bcast<int>(sizes[n], static_cast<int>(n));
    }

    // Now that we know how many items each node has, we can determine the
    // displacements
    std::vector<int> disps(static_cast<std::size_t>(size), 0);
    int disps_counter = 0;
    for (std::size_t n = 0; n < static_cast<std::size_t>(size); n++) {
      disps[n] = disps_counter;
      disps_counter += sizes[n];
    }

    // Temporary vector to recieve result on the root
    std::vector<T> tmp_rcv;
    if (rank == root) {
      std::size_t Ntot = 0;
      for (std::size_t n = 0; n < static_cast<std::size_t>(size); n++) {
        Ntot += static_cast<std::size_t>(sizes[n]);
      }
      tmp_rcv.resize(Ntot);
    }

    // Do the Gatherv
    int err =
        MPI_Gatherv(&vals[0], static_cast<int>(vals.size()), dtype<T>(),
                    &tmp_rcv[0], &sizes[0], &disps[0], dtype<T>(), root, com);
    check_error(err, file, line);

    if (rank == root) {
      vals = tmp_rcv;
    } else {
      vals.clear();
    }
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)file;
  (void)line;
#endif
}

template <typename T>
void Scatterv(std::vector<T>& vals, int root, const char* file = __FILE__,
              int line = __LINE__) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    // First, we need to know how many values there are on root to send,
    // and then how many values each node should get
    uint64_t Ntot = vals.size();
    Bcast<uint64_t>(Ntot, root);

    const int base = static_cast<int>(Ntot) / size;
    const int remainder = static_cast<int>(Ntot) - (size * base);

    std::vector<int> sizes(static_cast<std::size_t>(size), 0);
    for (std::size_t n = 0; n < static_cast<std::size_t>(size); n++) {
      sizes[n] = base;
      if (static_cast<int>(n) < remainder) sizes[n]++;
    }

    // Now we calculate the displacements if we are root
    std::vector<int> disps;
    if (rank == root) {
      disps.resize(static_cast<std::size_t>(size), 0);

      int disps_counter = 0;
      for (std::size_t n = 0; n < static_cast<std::size_t>(size); n++) {
        disps[n] = disps_counter;
        disps_counter += sizes[n];
      }
    }

    // Make a receiving buffer
    std::vector<T> tmp_rcv;
    tmp_rcv.resize(
        static_cast<std::size_t>(sizes[static_cast<std::size_t>(rank)]));

    // Do the Scatterv
    int err = MPI_Scatterv(&vals[0], &sizes[0], &disps[0], dtype<T>(),
                           &tmp_rcv[0], sizes[static_cast<std::size_t>(rank)],
                           dtype<T>(), root, com);
    check_error(err, file, line);

    vals = tmp_rcv;
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)file;
  (void)line;
#endif
}

}  // namespace mpi

#endif
