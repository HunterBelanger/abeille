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
#include <simulation/particle.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>

#include <type_traits>
#include <utility>

#ifdef ABEILLE_USE_MPI
#include <mpi.h>
#endif

namespace mpi {

Timer timer;

#ifdef ABEILLE_USE_MPI
const Com com = MPI_COMM_WORLD;
const DType Bool = MPI_C_BOOL;
const DType Int = MPI_INT;
const DType Double = MPI_DOUBLE;
const DType UInt16 = MPI_UINT16_T;
const DType UInt32 = MPI_UINT32_T;
const DType UInt64 = MPI_UINT64_T;
DType BParticle;
DType KeyType;
DType KeyUInt32Pair;

const OpType Sum = MPI_SUM;
const OpType And = MPI_LAND;
const OpType Or = MPI_LOR;
#endif

std::vector<uint64_t> node_nparticles;
std::vector<uint64_t> node_nparticles_noise;

int size = 1;
int rank = 0;

void register_banked_particle_type() {
#ifdef ABEILLE_USE_MPI
  // Ensure that BankedParticle is standard layout ! This is
  // required by MPI.
  static_assert(std::is_standard_layout<BankedParticle>::value);

  BankedParticle p;
  int err = 0;
  constexpr std::size_t BP_NUM_MEMBERS = 14;
  int sizes[BP_NUM_MEMBERS]{3, 3, 1, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1};
  DType dtypes[BP_NUM_MEMBERS]{Double, Double, Double, Double, Double,
                               UInt64, UInt64, UInt64, Bool,   Double,
                               Double, Double, Double, Double};
  MPI_Aint disps[BP_NUM_MEMBERS];
  MPI_Get_address(&p.r, &disps[0]);
  MPI_Get_address(&p.u, &disps[1]);
  MPI_Get_address(&p.E, &disps[2]);
  MPI_Get_address(&p.wgt, &disps[3]);
  MPI_Get_address(&p.wgt2, &disps[4]);
  MPI_Get_address(&p.parent_history_id, &disps[5]);
  MPI_Get_address(&p.parent_daughter_id, &disps[6]);
  MPI_Get_address(&p.family_id, &disps[7]);
  MPI_Get_address(&p.parents_previous_was_virtual, &disps[8]);
  MPI_Get_address(&p.parents_previous_position, &disps[9]);
  MPI_Get_address(&p.parents_previous_direction, &disps[10]);
  MPI_Get_address(&p.parents_previous_previous_energy, &disps[11]);
  MPI_Get_address(&p.parents_previous_energy, &disps[12]);
  MPI_Get_address(&p.Esmp_parent, &disps[13]);
  for (int i = BP_NUM_MEMBERS - 1; i >= 0; --i) {
    disps[i] -= disps[0];
  }

  err =
      MPI_Type_create_struct(BP_NUM_MEMBERS, sizes, disps, dtypes, &BParticle);
  check_error(err);

  err = MPI_Type_commit(&BParticle);
  check_error(err);
#endif
}

void initialize_mpi(int* argc, char*** argv) {
  timer.reset();

#ifdef ABEILLE_USE_MPI
  int err = 0;

  // Initialize MPI
  err = MPI_Init(argc, argv);
  check_error(err);

  // Get worls size and our rank
  err = MPI_Comm_size(com, &size);
  check_error(err);

  err = MPI_Comm_rank(com, &rank);
  check_error(err);

  // Register required custom types.
  // !!! MUST BE DONE IN CERTAIN ORDER !!!
  register_banked_particle_type();
#else
  (void)argc;
  (void)argv;
#endif
}

void finalize_mpi() {
#ifdef ABEILLE_USE_MPI
  MPI_Finalize();
#endif
}

void abort_mpi() {
#ifdef ABEILLE_USE_MPI
  MPI_Abort(com, 1);
#endif
}

void synchronize() {
#ifdef ABEILLE_USE_MPI
  if (size > 1) {
    timer.start();
    MPI_Barrier(com);
    timer.stop();
  }
#endif
}

void check_error(int err) {
#ifdef ABEILLE_USE_MPI
  if (err != MPI_SUCCESS) {
    // First, we should get the error string
    char err_str[MPI_MAX_ERROR_STRING];
    int str_len = 0;
    MPI_Error_string(err, err_str, &str_len);
    std::string mssg(err_str, static_cast<std::size_t>(str_len));
    fatal_error(mssg);
  }
#else
  (void)err;
#endif
}
}  // namespace mpi