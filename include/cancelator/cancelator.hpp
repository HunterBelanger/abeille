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
#ifndef CANCELATOR_H
#define CANCELATOR_H

#include <simulation/particle.hpp>
#include <simulation/particle_mover.hpp>
#include <utils/rng.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>

class Cancelator {
 public:
  virtual ~Cancelator() = default;

  virtual bool add_particle(BankedParticle& p) = 0;
  virtual void perform_cancellation() = 0;
  virtual std::vector<BankedParticle> get_new_particles(RNG& rng) = 0;
  virtual void clear() = 0;
  virtual void check_particle_mover_compatibility(
      const std::shared_ptr<IParticleMover>& pmover) const = 0;

  bool cancel_dual_weights() const { return dual_weights_; }
  void set_cancel_dual_weights(bool dw) { dual_weights_ = dw; }

 protected:
  bool dual_weights_ {false};
};

std::shared_ptr<Cancelator> make_cancelator(const YAML::Node& node);

#endif
