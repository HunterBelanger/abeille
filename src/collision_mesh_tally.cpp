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
#include <simulation/collision_mesh_tally.hpp>
#include <sstream>
#include <utils/error.hpp>
#include <utils/position.hpp>
#include <utils/settings.hpp>
#include <vector>

void CollisionMeshTally::score_collision(const Particle& p,
                                         MaterialHelper& mat) {
  double Et = mat.Et(p.E());

  double scr = 1. / (Et * net_weight);

  /*if (p.wgt() == 0. && p.wgt2() == 0.) {
    std::stringstream out;
    out << " Zero collision weights.\n";
    std::cout << out.str();
  }*/

  int i = static_cast<int>(std::floor((p.r().x() - r_low.x()) / dx));
  int j = static_cast<int>(std::floor((p.r().y() - r_low.y()) / dy));
  int k = static_cast<int>(std::floor((p.r().z() - r_low.z()) / dz));
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
        scr *= p.wgt() * mat.Emt(mt, p.E());
        break;

      case Quantity::RealFlux:
        scr *= p.wgt();
        break;

      case Quantity::ImgFlux:
        scr *= p.wgt2();
        break;
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    tally_gen(uE, ui, uj, uk) += scr;
  }
}

std::string CollisionMeshTally::quantity_str() const {
  switch (quantity) {
    case Quantity::Flux:
      return "Flux";
      break;

    case Quantity::Elastic:
      return "Elastic";
      break;

    case Quantity::Absorption:
      return "Absorption";
      break;

    case Quantity::Fission:
      return "Fission\n";
      break;

    case Quantity::Total:
      return "Total";
      break;

    case Quantity::MT:
      return "MT = " + std::to_string(mt);
      break;

    case Quantity::RealFlux:
      return "RealFlux";
      break;

    case Quantity::ImgFlux:
      return "ImgFlux";
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
    std::string mssg = "Now valid low entry for mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  xl = node["low"][0].as<double>();
  yl = node["low"][1].as<double>();
  zl = node["low"][2].as<double>();
  Position plow(xl, yl, zl);

  // Get hi position
  double xh = 0., yh = 0., zh = 0.;
  if (!node["hi"] || !node["hi"].IsSequence() || node["hi"].size() != 3) {
    std::string mssg = "Now valid hi entry for mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  xh = node["hi"][0].as<double>();
  yh = node["hi"][1].as<double>();
  zh = node["hi"][2].as<double>();
  Position phi(xh, yh, zh);

  // Check positions
  if (plow.x() >= phi.x() || plow.y() >= phi.y() || plow.x() >= phi.x()) {
    std::string mssg =
        "Low coordinates for mesh tally must be less than hi coordinates.";
    fatal_error(mssg, __FILE__, __LINE__);
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
      std::string mssg =
          "Mesh tally shapes must be values greater than or equal to 1.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    nx = static_cast<uint64_t>(tmp_nx);
    ny = static_cast<uint64_t>(tmp_ny);
    nz = static_cast<uint64_t>(tmp_nz);
  }

  // Get energy bounds
  std::vector<double> ebounds;
  if (!node["energy-bounds"] || !node["energy-bounds"].IsSequence()) {
    std::string mssg = "No valid energy-bounds enetry provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  ebounds = node["energy-bounds"].as<std::vector<double>>();
  if (ebounds.size() < 2) {
    std::string mssg = "Energy-bounds must have at least two entries.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else if (!std::is_sorted(ebounds.begin(), ebounds.end())) {
    std::string mssg = "Energy-bounds must be sorted.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else if (ebounds.front() < 0.) {
    std::string mssg = "All energy-bounds entries must be positive.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get name
  std::string fname;
  if (!node["name"] || !node["name"].IsScalar()) {
    std::string mssg = "No valid name provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  fname = node["name"].as<std::string>();

  // Get tally quantity
  uint32_t mt = 0;
  std::string quant_str = "none";
  Quantity quantity = Quantity::Flux;
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    std::string mssg = "No quantity entry provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
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
      std::string mssg = "Cannot do MT tallies in multi-group mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Check for mt
    if (!node["mt"] || !node["mt"].IsScalar()) {
      std::string mssg =
          "Quantity of \"mt\" selected, but no provided mt value.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    int32_t tmp_mt = node["mt"].as<int32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      std::string mssg =
          "The value " + std::to_string(tmp_mt) + " is not a valid MT.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "real-flux") {
    quantity = Quantity::RealFlux;
  } else if (quant_str == "img-flux") {
    quantity = Quantity::ImgFlux;
  } else {
    std::string mssg = "Unkown tally quantity \"" + quant_str + "\".";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Construct based on estimator type
  return std::make_shared<CollisionMeshTally>(plow, phi, nx, ny, nz, ebounds,
                                              quantity, fname, mt);
}
