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
#include <simulation/collision_operators/branching_collision.hpp>
#include <simulation/collision_operators/branchless_isotope_collision.hpp>
#include <simulation/collision_operators/branchless_material_collision.hpp>
#include <simulation/collision_operators/fission_bank_saver.hpp>
#include <simulation/collision_operators/fixed_source_fission_operator.hpp>
#include <simulation/collision_operators/keff_fission_operator.hpp>
#include <simulation/collision_operators/noise_branching_collision.hpp>
#include <simulation/collision_operators/secondary_bank_saver.hpp>
#include <simulation/particle_mover.hpp>
#include <simulation/transport_operators/carter_tracker.hpp>
#include <simulation/transport_operators/delta_tracker.hpp>
#include <simulation/transport_operators/implicit_leakage_delta_tracker.hpp>
#include <simulation/transport_operators/surface_tracker.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <string>
#include <variant>

using TransporterObjct =
    std::variant<SurfaceTracker, DeltaTracker, ImplicitLeakageDeltaTracker,
                 CarterTracker>;

TransporterObjct make_transporter(const YAML::Node& sim) {
  // Get the transporter
  if (!sim["transport-operator"]) {
    // No transport operator entry, use SurfaceTracking by default
    return SurfaceTracker();
  } else if (sim["transport-operator"].IsMap() == false) {
    fatal_error("Invalid transport-operator entry in simulation.");
  }

  const YAML::Node& to = sim["transport-operator"];

  // Get the transporter type
  if (!to["type"] || to["type"].IsScalar() == false) {
    fatal_error("No type entry in transport-operator");
  }
  std::string type = to["type"].as<std::string>();

  // Check for the sampling xs data
  std::vector<double> xs, egrid;
  if (type == "carter-tracking") {
    if (!to["sampling-xs"] || to["sampling-xs"].IsMap() == false) {
      fatal_error("No sampling-xs entry present for cater-tracker.");
    }
    const YAML::Node& sxs = to["sampling-xs"];

    // Get the energy grid for CE
    if (settings::energy_mode == settings::EnergyMode::CE) {
      if (!sxs["energy"] || sxs["energy"].IsSequence() == false) {
        fatal_error("No energy entry in sampling-xs.");
      }
      egrid = sxs["energy"].as<std::vector<double>>();
    }

    // Get the xs values
    if (!sxs["xs"] || sxs["xs"].IsSequence() == false) {
      fatal_error("No xs entry in sampling-xs.");
    }
    xs = sxs["xs"].as<std::vector<double>>();
  }

  if (type == "surface-tracking") {
    return SurfaceTracker();
  } else if (type == "delta-tracking") {
    return DeltaTracker();
  } else if (type == "implicit-leakage-delta-tracking") {
    return ImplicitLeakageDeltaTracker();
  } else if (type == "carter-tracking") {
    if (settings::energy_mode == settings::EnergyMode::CE) {
      return CarterTracker(egrid, xs);
    } else {
      return CarterTracker(xs);
    }
  } else {
    fatal_error("Unknown transport-operator type " + type + ".");
  }

  // NEVER GETS HERE
  return SurfaceTracker();
}

template <class CollOp>
std::shared_ptr<IParticleMover> make_concrete_mover(const YAML::Node& sim,
                                                    const CollOp& cop) {
  TransporterObjct top = make_transporter(sim);

  if (std::holds_alternative<SurfaceTracker>(top)) {
    return std::make_shared<ParticleMover<SurfaceTracker, CollOp>>(
        std::get<SurfaceTracker>(top), cop);
  } else if (std::holds_alternative<DeltaTracker>(top)) {
    return std::make_shared<ParticleMover<DeltaTracker, CollOp>>(
        std::get<DeltaTracker>(top), cop);
  } else if (std::holds_alternative<ImplicitLeakageDeltaTracker>(top)) {
    return std::make_shared<ParticleMover<ImplicitLeakageDeltaTracker, CollOp>>(
        std::get<ImplicitLeakageDeltaTracker>(top), cop);
  } else if (std::holds_alternative<CarterTracker>(top)) {
    return std::make_shared<ParticleMover<CarterTracker, CollOp>>(
        std::get<CarterTracker>(top), cop);
  }

  // Should never get here !
  fatal_error("Something has gone very wrong");
  return nullptr;
}

std::shared_ptr<INoiseParticleMover> make_concrete_noise_mover(const YAML::Node& sim, const NoiseBranchingCollision& cop) {
  TransporterObjct top = make_transporter(sim);

  if (std::holds_alternative<SurfaceTracker>(top)) {
    return std::make_shared<NoiseParticleMover<SurfaceTracker>>(std::get<SurfaceTracker>(top), cop);
  } else if (std::holds_alternative<DeltaTracker>(top)) {
    return std::make_shared<NoiseParticleMover<DeltaTracker>>(std::get<DeltaTracker>(top), cop);
  } else if (std::holds_alternative<ImplicitLeakageDeltaTracker>(top)) {
    return std::make_shared<NoiseParticleMover<ImplicitLeakageDeltaTracker>>(std::get<ImplicitLeakageDeltaTracker>(top), cop);
  } else if (std::holds_alternative<CarterTracker>(top)) {
    return std::make_shared<NoiseParticleMover<CarterTracker>>(std::get<CarterTracker>(top), cop);
  }

  // Should never get here !
  fatal_error("Something has gone very wrong");
  return nullptr;
}

std::shared_ptr<INoiseParticleMover> make_noise_particle_mover(const YAML::Node& sim) {
  bool noise_generations = true;
  if (sim["noise-generations"] && sim["noise-generations"].IsScalar()) {
    noise_generations = sim["noise-generations"].as<bool>();
  } else if (sim["noise-generations"]) {
    fatal_error("Invalid noise-generations entry in noise simulation.");
  }
  
  // If this is for a noise simulation, we only have one type of collison
  // operator at the moment, so nothing to check other than the mode.
  return make_concrete_noise_mover(sim, NoiseBranchingCollision(noise_generations));
}

std::shared_ptr<IParticleMover> make_particle_mover(const YAML::Node& sim,
                                                    settings::SimMode mode) {
  if (mode == settings::SimMode::NOISE) {
    fatal_error("Cannot create a particle mover for noise here."); 
  }

  // Next, we get the collision operator type
  if (!sim["collision-operator"]) {
    // No collision operator entry, use branching by default
    if (mode == settings::SimMode::FIXED_SOURCE) {
      auto cop = BranchingCollision<FixedSourceFissionOperator>();
      return make_concrete_mover(sim, cop);
    } else {
      auto cop = BranchingCollision<KeffFissionOperator>();
      return make_concrete_mover(sim, cop);
    }
  } else if (sim["collision-operator"].IsMap() == false) {
    fatal_error("Invalid collision-operator entry in simulation.");
  }

  const YAML::Node& co = sim["collision-operator"];
  if (!co["type"] || co["type"].IsScalar() == false) {
    fatal_error("No type entry in collision-operator");
  }
  std::string type = co["type"].as<std::string>();

  // Check for a splitting entry (present with branchless)
  bool splitting = false;
  if (co["splitting"] && co["splitting"].IsScalar()) {
    splitting = co["splitting"].as<bool>();
  } else if (co["splitting"]) {
    fatal_error("Invalid splitting entry in collision-operator.");
  }

  if (splitting && type == "branching") {
    fatal_error("Cannot use splitting with branchless collisions.");
  }

  if (type == "branching") {
    if (mode == settings::SimMode::FIXED_SOURCE) {
      auto cop = BranchingCollision<FixedSourceFissionOperator>();
      return make_concrete_mover(sim, cop);
    } else {
      auto cop = BranchingCollision<KeffFissionOperator>();
      return make_concrete_mover(sim, cop);
    }
  } else if (type == "branchless-isotope") {
    if (mode == settings::SimMode::FIXED_SOURCE) {
      auto cop = BranchlessIsotopeCollision<SecondaryBankSaver>();
      cop.set_splitting(splitting);
      return make_concrete_mover(sim, cop);
    } else {
      auto cop = BranchlessIsotopeCollision<FissionBankSaver>();
      cop.set_splitting(splitting);
      return make_concrete_mover(sim, cop);
    }
  } else if (type == "branchless-material") {
    if (mode == settings::SimMode::FIXED_SOURCE) {
      auto cop = BranchlessMaterialCollision<SecondaryBankSaver>();
      cop.set_splitting(splitting);
      return make_concrete_mover(sim, cop);
    } else {
      auto cop = BranchlessMaterialCollision<FissionBankSaver>();
      cop.set_splitting(splitting);
      return make_concrete_mover(sim, cop);
    }
  }

  // Should never get here !
  fatal_error("Quelque chose ne va pas de tout...");
  return nullptr;
}