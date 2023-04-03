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
#include <algorithm>
#include <simulation/track_length_mesh_tally.hpp>
#include <sstream>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/position.hpp>

inline double TrackLengthMeshTally::get_base_score(const Particle& p,
                                                   MaterialHelper& mat) const {
  // Calculate base score, absed on the quantity
  double base_score = 1. / net_weight;
  // Must multiply score by correct factor depending on the observed
  // quantity. q = 0 corresponds to just the flux, so no modification.
  switch (quantity) {
    case Quantity::Flux:
      base_score *= p.wgt();
      break;

    case Quantity::Elastic:
      base_score *= p.wgt() * mat.Eelastic(p.E());
      break;

    case Quantity::Absorption:
      base_score *= p.wgt() * mat.Ea(p.E());
      break;

    case Quantity::Fission:
      base_score *= p.wgt() * mat.Ef(p.E());
      break;

    case Quantity::Total:
      base_score *= p.wgt() * mat.Et(p.E());
      break;

    case Quantity::MT:
      base_score *= p.wgt() * mat.Emt(mt, p.E());
      break;

    case Quantity::RealFlux:
      base_score *= p.wgt();
      break;

    case Quantity::ImgFlux:
      base_score *= p.wgt2();
      break;
  }

  return base_score;
}

void TrackLengthMeshTally::score_flight(const Particle& p, double d,
                                        MaterialHelper& mat) {
  // Local position and direction copies which we will use for tallying
  Position r = p.r();
  Direction u = p.u();

  // Initialize the spatial indices with the starting position
  // of the particle
  int i = 0, j = 0, k = 0;
  std::array<int, 3> on{0, 0, 0};
  initialize_indices(r, u, i, j, k, on);

  // First, check if we are inside the tally region
  bool inside_tally_region = false;
  if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
      j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
    inside_tally_region = true;
  }

  // If we aren't inside, are we going to hit it along the trajectory ?
  // If so, we can fast-forward to that position.
  if (!inside_tally_region) {
    if (!find_entry_point(r, u, d)) {
      return;
    }

    // We should now be on the boundary of the mesh, so we need to
    // re-initialize our indicies, and check to make sure we are inside
    // the tally region.
    initialize_indices(r, u, i, j, k, on);
    if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
        j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
      inside_tally_region = true;
    } else {
      // This is a problem, in theory, we should now be inside the tally
      // region. We will therefore spew a warning here.
      warning("Could not locate tile after fast forward to mesh entry.\n");
    }
  }

  // Calculate base score, absed on the quantity
  double base_score = this->get_base_score(p, mat);

  // Get energy index with linear search
  int l = -1;
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
  uint64_t uE = static_cast<uint64_t>(l);

  // Distance remaining to tally
  double distance_remaining = d;

  while (distance_remaining > 0.) {
    // Distance we will travel in this cell
    auto next_tile = distance_to_next_index(r, u, i, j, k, on);

    if (next_tile.first == INF) {
      // Something went wrong.... Don't score.
      Output::instance()->save_warning("Problem encountered with mesh tally " +
                                       fname + ".");
      break;
    } else if (next_tile.first < 0.) {
      // Something went wrong.... Don't score.
      warning("Negative distance encountered with mesh tally.");
    }

    double d_tile = std::min(next_tile.first, distance_remaining);

    // Make the score if we are in a valid cell
    if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
        j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
      uint64_t ui = static_cast<uint64_t>(i);
      uint64_t uj = static_cast<uint64_t>(j);
      uint64_t uk = static_cast<uint64_t>(k);
#ifdef _OPENMP
#pragma omp atomic
#endif
      tally_gen(uE, ui, uj, uk) += d_tile * base_score;
    } else {
      // If we arrive here, it means that we have left the tally region
      // when were we initially inside it. We can return here, as it's
      // impossible to go back in.
      return;
    }

    // Remove the traveled distance
    distance_remaining -= d_tile;

    if (distance_remaining <= 0.) break;

    // Update the position and cell indices
    r = r + d_tile * u;
    update_indices(next_tile.second, i, j, k, on);

  }  // While we still have to travel
}

bool TrackLengthMeshTally::find_entry_point(Position& r, const Direction& u,
                                            double& d_flight) const {
  double d_min = (r_low.x() - r.x()) / u.x();
  double d_max = (r_hi.x() - r.x()) / u.x();

  if (d_min > d_max) {
    std::swap(d_min, d_max);
  }

  double d_y_min = (r_low.y() - r.y()) / u.y();
  double d_y_max = (r_hi.y() - r.y()) / u.y();

  if (d_y_min > d_y_max) {
    std::swap(d_y_min, d_y_max);
  }

  if ((d_min > d_y_max) || (d_y_min > d_max)) {
    return false;
  }

  if (d_y_min > d_min) {
    d_min = d_y_min;
  }

  if (d_y_max < d_max) {
    d_max = d_y_max;
  }

  double d_z_min = (r_low.z() - r.z()) / u.z();
  double d_z_max = (r_hi.z() - r.z()) / u.z();

  if (d_z_min > d_z_max) {
    std::swap(d_z_min, d_z_max);
  }

  if ((d_min > d_z_max) || (d_z_min > d_max)) {
    return false;
  }

  if (d_z_min > d_min) {
    d_min = d_z_min;
  }

  if (d_z_max < d_max) {
    d_max = d_z_max;
  }

  if (d_max < d_min) {
    std::swap(d_max, d_min);
  }

  if ((d_max < 0.) && (d_min < 0.)) {
    return false;
  }

  if (d_min < 0.) {
    // If we are here, this means that r is actually inside the mesh, but is
    // really close to the edge, and we have a direction taking us out.
    // We should return false here, so that we don't score anything for this
    // particle track.
    return false;
  }

  // If we get here, we intersect the box. Let's update the position and the
  // flight distance.
  r = r + d_min * u;
  d_flight -= d_min;
  return true;
}

void TrackLengthMeshTally::initialize_indices(const Position& r,
                                              const Direction& u, int& i,
                                              int& j, int& k,
                                              std::array<int, 3>& on) {
  i = static_cast<int>(std::floor((r.x() - r_low.x()) / dx));
  j = static_cast<int>(std::floor((r.y() - r_low.y()) / dy));
  k = static_cast<int>(std::floor((r.z() - r_low.z()) / dz));
  on.fill(0);

  // Must handle case of being on a tile boundary
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get tile boundaries
  double xl = xc - 0.5 * dx;
  double xh = xc + 0.5 * dx;
  double yl = yc - 0.5 * dy;
  double yh = yc + 0.5 * dy;
  double zl = zc - 0.5 * dz;
  double zh = zc + 0.5 * dz;

  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      i--;
      on[0] = 1;
    } else {
      on[0] = -1;
    }
  } else if (std::abs(xh - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      on[0] = 1;
    } else {
      i++;
      on[0] = -1;
    }
  }

  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      j--;
      on[1] = 1;
    } else {
      on[1] = -1;
    }
  } else if (std::abs(yh - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      on[1] = 1;
    } else {
      j++;
      on[1] = -1;
    }
  }

  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      k--;
      on[2] = 1;
    } else {
      on[2] = -1;
    }
  } else if (std::abs(zh - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      on[2] = 1;
    } else {
      k++;
      on[2] = -1;
    }
  }
}

std::pair<double, int> TrackLengthMeshTally::distance_to_next_index(
    const Position& r, const Direction& u, int i, int j, int k,
    const std::array<int, 3>& on) {
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get relative position in cell
  Position r_tile(r.x() - xc, r.y() - yc, r.z() - zc);

  // Set our initial value for the distance and the index change
  double dist = INF;
  int key = 0;

  // Check all six sides
  double diff_xl = -dx * 0.5 - r_tile.x();
  double diff_xh = dx * 0.5 - r_tile.x();
  double diff_yl = -dy * 0.5 - r_tile.y();
  double diff_yh = dy * 0.5 - r_tile.y();
  double diff_zl = -dz * 0.5 - r_tile.z();
  double diff_zh = dz * 0.5 - r_tile.z();

  double d_xl = diff_xl / u.x();
  double d_xh = diff_xh / u.x();
  double d_yl = diff_yl / u.y();
  double d_yh = diff_yh / u.y();
  double d_zl = diff_zl / u.z();
  double d_zh = diff_zh / u.z();

  if (d_xl > 0. && d_xl < dist && on[0] != -1) {
    dist = d_xl;
    key = -1;
  }

  if (d_xh > 0. && d_xh < dist && on[0] != 1) {
    dist = d_xh;
    key = 1;
  }

  if (d_yl > 0. && d_yl < dist && on[1] != -1) {
    dist = d_yl;
    key = -2;
  }

  if (d_yh > 0. && d_yh < dist && on[1] != 1) {
    dist = d_yh;
    key = 2;
  }

  if (d_zl > 0. && d_zl < dist && on[2] != -1) {
    dist = d_zl;
    key = -3;
  }

  if (d_zh > 0. && d_zh < dist && on[2] != 1) {
    dist = d_zh;
    key = 3;
  }

  return {dist, key};
}

void TrackLengthMeshTally::update_indices(int key, int& i, int& j, int& k,
                                          std::array<int, 3>& on) {
  // Must initially fill with zero, so that we don't stay on top
  // of other surfaces the entire time
  on.fill(0);

  switch (key) {
    case -1:
      i--;
      on[0] = 1;
      break;

    case 1:
      i++;
      on[0] = -1;
      break;

    case -2:
      j--;
      on[1] = 1;
      break;

    case 2:
      j++;
      on[1] = -1;
      break;

    case -3:
      k--;
      on[2] = 1;
      break;

    case 3:
      k++;
      on[2] = -1;
      break;

    default:
      break;
  }
}

std::string TrackLengthMeshTally::quantity_str() const {
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

  // Never gets here
  return "unkown";
}

std::shared_ptr<TrackLengthMeshTally> make_track_length_mesh_tally(
    const YAML::Node& node) {
  using Quantity = TrackLengthMeshTally::Quantity;

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
      std::stringstream mssg;
      mssg << "The value " << tmp_mt << " is not a valid MT.";
      fatal_error(mssg.str());
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "real-flux") {
    quantity = Quantity::RealFlux;
  } else if (quant_str == "img-flux") {
    quantity = Quantity::ImgFlux;
  } else {
    fatal_error("Unkown tally quantity \"" + quant_str + "\".");
  }

  if ((settings::tracking == settings::TrackingMode::DELTA_TRACKING ||
       settings::tracking == settings::TrackingMode::CARTER_TRACKING) &&
      (quantity != Quantity::Flux && quantity != Quantity::RealFlux &&
       quantity != Quantity::ImgFlux)) {
    fatal_error(
        "Cannot use track-length estimators for non-flux quantities with "
        "delta-tracking or carter-tracking.");
  }

  return std::make_shared<TrackLengthMeshTally>(plow, phi, nx, ny, nz, ebounds,
                                                quantity, fname, mt);
}
