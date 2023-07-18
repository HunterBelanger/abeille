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
#ifndef MESH_TALLY_H
#define MESH_TALLY_H

#include <materials/material_helper.hpp>
#include <materials/nuclide.hpp>
#include <simulation/particle.hpp>

#include <yaml-cpp/yaml.h>
#include <ndarray.hpp>

#include <memory>

class MeshTally {
 public:
  enum class Quantity {
    Flux,
    Total,
    Elastic,
    Absorption,
    Fission,
    MT,
    RealFlux,
    ImgFlux
  };

  MeshTally(Position low, Position hi, uint64_t nx, uint64_t ny, uint64_t nz,
            const std::vector<double>& ebounds, std::string fname);
  virtual ~MeshTally() = default;

  virtual std::string estimator_str() const = 0;

  virtual std::string quantity_str() const = 0;

  virtual std::uint32_t mt() const = 0;

  void set_net_weight(double W);

  // multiplier must be the same across ALL MPI processes, as it is
  // applied after the MPI reduction of the generation scores.
  void record_generation(double multiplier = 1.);

  void clear_generation();

  void write_tally();

 protected:
  Position r_low, r_hi;
  uint64_t Nx, Ny, Nz, g;
  double dx, dy, dz, dx_inv, dy_inv, dz_inv, net_weight;
  std::vector<double> energy_bounds;
  std::string fname;

  NDArray<double> tally_gen;
  NDArray<double> tally_avg;
  NDArray<double> tally_var;
};

#endif
