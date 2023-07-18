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
#ifndef SURFACE_TRACKER_H
#define SURFACE_TRACKER_H

#include <materials/nuclide.hpp>
#include <simulation/transporter.hpp>

class SurfaceTracker : public Transporter {
 public:
  SurfaceTracker(std::shared_ptr<Tallies> i_t) : Transporter{i_t} {}
  ~SurfaceTracker() = default;

  std::vector<BankedParticle> transport(
      std::vector<Particle>& bank, bool noise = false,
      std::vector<BankedParticle>* noise_bank = nullptr,
      const NoiseMaker* noise_maker = nullptr);

};  // SurfaceTracker

#endif  // MG_SURFACE_TRACKER_H
