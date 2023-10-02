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

#include <simulation/cancelator.hpp>

#include <unordered_map>

class ApproximateMeshCancelator : public Cancelator {
 public:
   std::unordered_map<int, std::vector<BankedParticle*>> bins; //MAKE PRIVATE LATER  

  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz);
  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz, std::vector<double> energy_bounds);

  void sync_keys(std::vector<int>& keys); // MAKE PRIVATE LATER
  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation(pcg32& rng) override final;
  void perform_cancellation_loop(pcg32& rng);
  void perform_cancellation_vector(pcg32& rng);
  std::vector<BankedParticle> get_new_particles(pcg32& rng) override final;
  void clear() override final;


 private:
  Position r_low, r_hi;
  std::vector<double> energy_edges;
  std::array<uint32_t, 4> shape;
  double dx, dy, dz;
};

std::shared_ptr<ApproximateMeshCancelator> make_approximate_mesh_cancelator(
    const YAML::Node& node);

#endif
