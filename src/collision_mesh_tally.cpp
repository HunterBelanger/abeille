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
#include <tallies/collision_mesh_tally.hpp>
#include <utils/error.hpp>
#include <utils/position.hpp>
#include <utils/settings.hpp>

#include <sstream>
#include <vector>

void CollisionMeshTally::score_collision(const Particle& p,
                                         MaterialHelper& mat) {
  double Et = mat.Et(p.E());

  double scr = 1. / (Et * net_weight);

  int i = static_cast<int>(std::floor((p.r().x() - r_low.x()) * dx_inv));
  int j = static_cast<int>(std::floor((p.r().y() - r_low.y()) * dy_inv));
  int k = static_cast<int>(std::floor((p.r().z() - r_low.z()) * dz_inv));
  int l = -1;

  // Get energy index with linear search
  for (size_t e = 0; e < energy_bounds.size() - 1; e++) {
    if (energy_bounds[e] <= p.E() && p.E() <= energy_bounds[e + 1]) {
      l = static_cast<int>(e);
      break;
    }
  }

  // Don't score anything if we didn't find an energy bin
  if (l == -1) {
    return;
  }

  if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
      j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
    uint64_t ui = static_cast<uint64_t>(i);
    uint64_t uj = static_cast<uint64_t>(j);
    uint64_t uk = static_cast<uint64_t>(k);
    uint64_t uE = static_cast<uint64_t>(l);

    // Must multiply score by correct factor depending on the observed
    // quantity. q = 0 corresponds to just the flux, so no modification.
    switch (quantity) {
      case Quantity::Flux:
        scr *= p.wgt();
        break;

      case Quantity::Elastic:
        scr *= p.wgt() * mat.Eelastic(p.E());
        break;

      case Quantity::Absorption:
        scr *= p.wgt() * mat.Ea(p.E());
        break;

      case Quantity::Fission:
        scr *= p.wgt() * mat.Ef(p.E());
        break;

      case Quantity::Total:
        scr *= p.wgt() * Et;
        break;

      case Quantity::MT:
        scr *= p.wgt() * mat.Emt(mt_, p.E());
        break;

      case Quantity::RealFlux:
        scr *= p.wgt();
        break;

      case Quantity::ImgFlux:
        scr *= p.wgt2();
        break;
    }
#ifdef ABEILLE_USE_OMP
#pragma omp atomic
#endif
    tally_gen(uE, ui, uj, uk) += scr;
  }
}

std::string CollisionMeshTally::quantity_str() const {
  switch (quantity) {
    case Quantity::Flux:
      return "flux";
      break;

    case Quantity::Elastic:
      return "elastic";
      break;

    case Quantity::Absorption:
      return "absorption";
      break;

    case Quantity::Fission:
      return "fission\n";
      break;

    case Quantity::Total:
      return "total";
      break;

    case Quantity::MT:
      return "mt";
      break;

    case Quantity::RealFlux:
      return "real-flux";
      break;

    case Quantity::ImgFlux:
      return "imag-flux";
      break;
  }

  // NEVER GETS HERE!
  return "unkown";
}

std::shared_ptr<CollisionMeshTally> make_collision_mesh_tally(
    const YAML::Node& node) {
  using Quantity = CollisionMeshTally::Quantity;

  // Get low position
  double xl = 0., yl = 0., zl = 0.;
  if (!node["low"] || !node["low"].IsSequence() || node["low"].size() != 3) {
    fatal_error("Now valid low entry for mesh tally.");
  }
  xl = node["low"][0].as<double>();
  yl = node["low"][1].as<double>();
  zl = node["low"][2].as<double>();
  Position plow(xl, yl, zl);

  // Get hi position
  double xh = 0., yh = 0., zh = 0.;
  if (!node["hi"] || !node["hi"].IsSequence() || node["hi"].size() != 3) {
    fatal_error("Now valid hi entry for mesh tally.");
  }
  xh = node["hi"][0].as<double>();
  yh = node["hi"][1].as<double>();
  zh = node["hi"][2].as<double>();
  Position phi(xh, yh, zh);

  // Check positions
  if (plow.x() >= phi.x() || plow.y() >= phi.y() || plow.x() >= phi.x()) {
    fatal_error(
        "Low coordinates for mesh tally must be less than hi coordinates.");
  }

  // Get shape
  uint64_t nx = 1, ny = 1, nz = 1;
  if (node["shape"] &&
      (!node["shape"].IsSequence() || node["shape"].size() != 3)) {
    std::string mssg = "Invalid shape provided to mesh tally.";
  } else if (node["shape"]) {
    int64_t tmp_nx = node["shape"][0].as<int64_t>();
    int64_t tmp_ny = node["shape"][1].as<int64_t>();
    int64_t tmp_nz = node["shape"][2].as<int64_t>();

    if (tmp_nx < 1 || tmp_ny < 1 || tmp_nz < 1) {
      fatal_error(
          "Mesh tally shapes must be values greater than or equal to 1.");
    }

    nx = static_cast<uint64_t>(tmp_nx);
    ny = static_cast<uint64_t>(tmp_ny);
    nz = static_cast<uint64_t>(tmp_nz);
  }

  // Get energy bounds
  std::vector<double> ebounds;
  if (!node["energy-bounds"] || !node["energy-bounds"].IsSequence()) {
    fatal_error("No valid energy-bounds enetry provided to mesh tally.");
  }
  ebounds = node["energy-bounds"].as<std::vector<double>>();
  if (ebounds.size() < 2) {
    fatal_error("Energy-bounds must have at least two entries.");
  } else if (!std::is_sorted(ebounds.begin(), ebounds.end())) {
    fatal_error("Energy-bounds must be sorted.");
  } else if (ebounds.front() < 0.) {
    fatal_error("All energy-bounds entries must be positive.");
  }

  // Get name
  std::string fname;
  if (!node["name"] || !node["name"].IsScalar()) {
    fatal_error("No valid name provided to mesh tally.");
  }
  fname = node["name"].as<std::string>();

  // Get tally quantity
  uint32_t mt = 0;
  std::string quant_str = "none";
  Quantity quantity = Quantity::Flux;
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    fatal_error("No quantity entry provided to mesh tally.");
  }
  quant_str = node["quantity"].as<std::string>();
  if (quant_str == "flux") {
    quantity = Quantity::Flux;
  } else if (quant_str == "total") {
    quantity = Quantity::Total;
  } else if (quant_str == "elastic") {
    quantity = Quantity::Elastic;
  } else if (quant_str == "absorption") {
    quantity = Quantity::Absorption;
  } else if (quant_str == "fission") {
    quantity = Quantity::Fission;
  } else if (quant_str == "mt") {
    quantity = Quantity::MT;

    if (settings::energy_mode == settings::EnergyMode::MG) {
      // Can't do an MT tally in MG mode !
      fatal_error("Cannot do MT tallies in multi-group mode.");
    }

    // Check for mt
    if (!node["mt"] || !node["mt"].IsScalar()) {
      fatal_error("Quantity of \"mt\" selected, but no provided mt value.");
    }

    int32_t tmp_mt = node["mt"].as<int32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      fatal_error("The value " + std::to_string(tmp_mt) +
                  " is not a valid MT.");
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "real-flux") {
    quantity = Quantity::RealFlux;
  } else if (quant_str == "imag-flux") {
    quantity = Quantity::ImgFlux;
  } else {
    fatal_error("Unkown tally quantity \"" + quant_str + "\".");
  }

  // Construct based on estimator type
  return std::make_shared<CollisionMeshTally>(plow, phi, nx, ny, nz, ebounds,
                                              quantity, fname, mt);
}
