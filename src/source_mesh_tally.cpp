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
#include <simulation/source_mesh_tally.hpp>
#include <utils/error.hpp>
#include <utils/position.hpp>
#include <utils/settings.hpp>

void SourceMeshTally::score_source(const BankedParticle& p) {
  int i = static_cast<int>(std::floor((p.r.x() - r_low.x()) / dx));
  int j = static_cast<int>(std::floor((p.r.y() - r_low.y()) / dy));
  int k = static_cast<int>(std::floor((p.r.z() - r_low.z()) / dz));
  int l = -1;

  // Get energy index with linear search
  for (size_t e = 0; e < energy_bounds.size() - 1; e++) {
    if (energy_bounds[e] <= p.E && p.E <= energy_bounds[e + 1]) {
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
    double scr = 1. / net_weight;
    switch (quantity) {
      case Quantity::Source:
        scr *= p.wgt;
        break;

      case Quantity::RealSource:
        scr *= p.wgt;
        break;

      case Quantity::ImagSource:
        scr *= p.wgt2;
        break;
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    tally_gen(uE, ui, uj, uk) += scr;
  }
}

std::shared_ptr<SourceMeshTally> make_source_mesh_tally(
    const YAML::Node& node) {
  using Quantity = SourceMeshTally::Quantity;

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
  std::string quant_str = "none";
  Quantity quantity = Quantity::Source;
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    fatal_error("No quantity entry provided to mesh tally.");
  }
  quant_str = node["quantity"].as<std::string>();
  if (quant_str == "source") {
    quantity = Quantity::Source;
  } else if (quant_str == "real-source") {
    quantity = Quantity::RealSource;
  } else if (quant_str == "img-source") {
    quantity = Quantity::ImagSource;
  } else {
    fatal_error("Unkown tally quantity \"" + quant_str + "\".");
  }

  return std::make_shared<SourceMeshTally>(plow, phi, nx, ny, nz, ebounds,
                                           quantity, fname);
}

std::string SourceMeshTally::quantity_str() const {
  switch (quantity) {
    case Quantity::Source:
      return "source";
      break;

    case Quantity::RealSource:
      return "real-source";
      break;

    case Quantity::ImagSource:
      return "img-source";
      break;
  }

  // Never get's here
  return "unknown";
}
