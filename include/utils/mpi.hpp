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
#ifndef MPI_H
#define MPI_H

#include <utils/error.hpp>
#include <utils/timer.hpp>

#include <cstdint>
#include <source_location>
#include <span>
#include <vector>

#ifdef ABEILLE_USE_MPI
#include <mpi.h>

#include <cancelator/exact_mg_cancelator.hpp>
#include <simulation/particle.hpp>
#endif

namespace mpi {

// Timer for all MPI operations
extern Timer timer;

#ifdef ABEILLE_USE_MPI
using Com = MPI_Comm;
using DType = MPI_Datatype;
using OpType = MPI_Op;

extern const Com com;
extern const DType Char;
extern const DType WChar;
extern const DType Bool;
extern const DType Int;
extern const DType Double;
extern const DType UInt16;
extern const DType UInt32;
extern const DType UInt64;
extern DType BParticle;
extern DType KeyType;
extern DType KeyUInt32Pair;

extern const OpType Sum;
extern const OpType And;
extern const OpType Or;

//==============================================================================
// MPI Data Types
template <typename T>
DType dtype();

template <>
inline DType dtype<char>() {
  return Char;
}

template <>
inline DType dtype<wchar_t>() {
  return WChar;
}

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
inline DType dtype<uint16_t>() {
  return UInt16;
}

template <>
inline DType dtype<uint32_t>() {
  return UInt32;
}

template <>
inline DType dtype<uint64_t>() {
  return UInt64;
}

template <>
inline DType dtype<BankedParticle>() {
  return BParticle;
}

template <>
inline DType dtype<ExactMGCancelator::Key>() {
  return KeyType;
}

template <>
inline DType dtype<std::pair<ExactMGCancelator::Key, uint32_t>>() {
  return KeyUInt32Pair;
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
void check_error(int err, const std::source_location& loc);

//==============================================================================
// MPI Operations
template <typename T>
void Bcast(T& val, int root,
           std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err = MPI_Bcast(&val, 1, dtype<T>(), root, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)val;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Bcast(std::vector<T>& vals, int root,
           std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();

    // First, we broadcast the size of the vector to all processes.
    std::size_t count = vals.size();
    Bcast(count, root);

    // Resize if necessary
    vals.resize(count);

    // Now we broadcast data.
    int err = MPI_Bcast(&vals[0], static_cast<int>(vals.size()), dtype<T>(),
                        root, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Send(T& val, int dest, int tag = 0,
          std::source_location loc = std::source_location::current()) {
  if (size < 2) {
    fatal_error("mpi::Send cannot be called with mpi::size less than 2", loc);
  }
#ifdef ABEILLE_USE_MPI
  timer.start();
  int err = MPI_Send(&val, 1, dtype<T>(), dest, tag, com);
  check_error(err, loc);
  timer.stop();
#else
  (void)val;
  (void)dest;
  (void)tag;
  (void)loc;
#endif
}

template <typename T>
void Send(std::span<T> vals, int dest, int tag = 0,
          std::source_location loc = std::source_location::current()) {
  if (size < 2) {
    fatal_error("mpi::Send cannot be called with mpi::size less than 2", loc);
  }
#ifdef ABEILLE_USE_MPI
  timer.start();
  int err = MPI_Send(&vals[0], static_cast<int>(vals.size()), dtype<T>(), dest,
                     tag, com);
  check_error(err, loc);
  timer.stop();
#else
  (void)vals;
  (void)dest;
  (void)tag;
  (void)loc;
#endif
}

template <typename T>
void Recv(T& val, int src, int tag = 0,
          std::source_location loc = std::source_location::current()) {
  if (size < 2) {
    fatal_error("mpi::Recv cannot be called with mpi::size less than 2", loc);
  }
#ifdef ABEILLE_USE_MPI
  timer.start();
  int err = MPI_Recv(&val, 1, dtype<T>(), src, tag, com, MPI_STATUS_IGNORE);
  check_error(err, loc);
  timer.stop();
#else
  (void)val;
  (void)src;
  (void)tag;
  (void)loc;
#endif
}

template <typename T>
void Recv(std::span<T> vals, int src, int tag = 0,
          std::source_location loc = std::source_location::current()) {
  if (size < 2) {
    fatal_error("mpi::Recv cannot be called with mpi::size less than 2", loc);
  }
#ifdef ABEILLE_USE_MPI
  timer.start();
  int err = MPI_Recv(&vals[0], static_cast<int>(vals.size()), dtype<T>(), src,
                     tag, com, MPI_STATUS_IGNORE);
  check_error(err, loc);
  timer.stop();
#else
  (void)vals;
  (void)src;
  (void)tag;
  (void)loc;
#endif
}

template <typename T>
void Allreduce_sum(T& val,
                   std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err = MPI_Allreduce(MPI_IN_PLACE, &val, 1, dtype<T>(), Sum, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)val;
  (void)loc;
#endif
}

template <typename T>
void Allreduce_sum(std::vector<T>& vals,
                   std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err =
        MPI_Allreduce(MPI_IN_PLACE, &vals[0], static_cast<int>(vals.size()),
                      dtype<T>(), Sum, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)vals;
  (void)loc;
#endif
}

template <typename T>
void Allreduce_sum(std::span<T> vals,
                   std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err =
        MPI_Allreduce(MPI_IN_PLACE, &vals[0], static_cast<int>(vals.size()),
                      dtype<T>(), Sum, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)vals;
  (void)loc;
#endif
}

template <typename T>
void Allreduce_or(T& val,
                  std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err = MPI_Allreduce(MPI_IN_PLACE, &val, 1, dtype<T>(), Or, com);
    check_error(err, loc);
    timer.stop();
  }
#else
  (void)val;
  (void)loc;
#endif
}

template <typename T>
void Reduce_sum(T& val, int root,
                std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    int err;
    if (mpi::rank == root) {
      err = MPI_Reduce(MPI_IN_PLACE, &val, 1, dtype<T>(), Sum, root, com);
    } else {
      err = MPI_Reduce(&val, nullptr, 1, dtype<T>(), Sum, root, com);
    }
    check_error(err, loc);

    timer.stop();
  }
#else
  (void)val;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Reduce_sum(std::vector<T>& vals, int root,
                std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();

    int err;
    if (mpi::rank == root) {
      err = MPI_Reduce(MPI_IN_PLACE, &vals[0], static_cast<int>(vals.size()),
                       dtype<T>(), Sum, root, com);
    } else {
      err = MPI_Reduce(&vals[0], nullptr, static_cast<int>(vals.size()),
                       dtype<T>(), Sum, root, com);
    }
    check_error(err, loc);

    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Reduce_sum(std::span<T> vals, int root,
                std::source_location loc = std::source_location::current()) {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();

    int err;
    if (mpi::rank == root) {
      err = MPI_Reduce(MPI_IN_PLACE, &vals[0], static_cast<int>(vals.size()),
                       dtype<T>(), Sum, root, com);
    } else {
      err = MPI_Reduce(&vals[0], nullptr, static_cast<int>(vals.size()),
                       dtype<T>(), Sum, root, com);
    }
    check_error(err, loc);

    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Gatherv(std::vector<T>& vals, int root,
             std::source_location loc = std::source_location::current()) {
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
    check_error(err, loc);

    if (rank == root) {
      vals.swap(tmp_rcv);
    } else {
      vals.clear();
    }
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)loc;
#endif
}

template <typename T>
void Scatterv(std::vector<T>& vals, int root,
              std::source_location loc = std::source_location::current()) {
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
    check_error(err, loc);

    vals = tmp_rcv;
    timer.stop();
  }
#else
  (void)vals;
  (void)root;
  (void)loc;
#endif
}

}  // namespace mpi

#endif
