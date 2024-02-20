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
#ifndef APPROXIMATE_MESH_CANCELATOR_H
#define APPROXIMATE_MESH_CANCELATOR_H

#include <cancelator/cancelator.hpp>

#include <unordered_map>

class ApproximateMeshCancelator : public Cancelator {
 public:
  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz, bool loop = false);
  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz, std::vector<double> energy_bounds,
                            bool loop = false);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation() override final;
  std::vector<BankedParticle> get_new_particles(RNG& rng) override final;
  void clear() override final;
  void check_particle_mover_compatibility(
      const std::shared_ptr<IParticleMover>& /*pmover*/) const override final {}

 private:
  std::unordered_map<int, std::vector<BankedParticle*>> bins;
  std::vector<double> energy_edges;
  std::array<uint32_t, 4> shape;
  Position r_low, r_hi;
  uint32_t Si, Sj, Sk, Sl;  // Strides for indexing
  double dx, dy, dz;
  bool loop;

  std::vector<int> sync_keys();
  void perform_cancellation_loop();
  void perform_cancellation_vector();
  void perform_cancellation_full_vector();
};

std::shared_ptr<ApproximateMeshCancelator> make_approximate_mesh_cancelator(
    const YAML::Node& node);

#endif
