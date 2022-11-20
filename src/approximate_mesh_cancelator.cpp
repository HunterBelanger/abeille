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
#include <cmath>
#include <simulation/approximate_mesh_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

ApproximateMeshCancelator::ApproximateMeshCancelator(Position low, Position hi,
                                                     uint32_t Nx, uint32_t Ny,
                                                     uint32_t Nz)
    : r_low(low),
      r_hi(hi),
      energy_edges(),
      shape{Nx, Ny, Nz, 1},
      dx(0.),
      dy(0.),
      dz(0.),
      bins() {
  // Make sure the low points are all lower than the high points
  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    std::string mssg =
        "Low position is not lower than hi position in "
        "ApproximateMeshCancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  dx = (r_hi.x() - r_low.x()) / static_cast<double>(shape[0]);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(shape[1]);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(shape[2]);
}

ApproximateMeshCancelator::ApproximateMeshCancelator(
    Position low, Position hi, uint32_t Nx, uint32_t Ny, uint32_t Nz,
    std::vector<double> energy_bounds)
    : r_low(low),
      r_hi(hi),
      energy_edges(energy_bounds),
      shape{Nx, Ny, Nz, 1},
      dx(0.),
      dy(0.),
      dz(0.),
      bins() {
  // Make sure the low points are all lower than the high points
  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    std::string mssg =
        "Low position is not lower than hi position in "
        "ApproximateMeshCancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  dx = (r_hi.x() - r_low.x()) / static_cast<double>(shape[0]);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(shape[1]);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(shape[2]);

  // Make sure energy bins are valid
  if (energy_edges.size() < 2) {
    std::string mssg =
        "energy_edges must have at least two entries in "
        "ApproximateMeshCancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(energy_edges.begin(), energy_edges.end())) {
    std::string mssg =
        "energy_edges must be sorted in ApproximateMeshCancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (energy_edges.front() < 0.) {
    std::string mssg =
        "All energy_edges must be greater than or equal to zero "
        "in ApproximateMeshCancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  shape[3] = static_cast<uint32_t>(energy_edges.size() - 1);
}

bool ApproximateMeshCancelator::add_particle(BankedParticle& p) {
  // Get bin indicies for spacial coordinates
  int i = std::floor((p.r.x() - r_low.x()) / dx);
  int j = std::floor((p.r.y() - r_low.y()) / dy);
  int k = std::floor((p.r.z() - r_low.z()) / dz);

  // Get energy index with linear search
  int l = -1;
  if (energy_edges.empty() == false) {
    for (size_t e = 0; e < energy_edges.size() - 1; e++) {
      if (energy_edges[e] <= p.E && p.E <= energy_edges[e + 1]) {
        l = static_cast<int>(e);
        break;
      }
    }
  } else {
    l = 0;
  }

  // If one of the indices is less than zero, don't keep the particle
  if (i < 0 || j < 0 || k < 0 || l < 0) {
    return false;
  }

  // if one of the indices is too large, don't keep the particle
  if (i >= static_cast<int>(shape[0]) || j >= static_cast<int>(shape[1]) ||
      k >= static_cast<int>(shape[2]) || l >= static_cast<int>(shape[3])) {
    return false;
  }

  // The particle fits into the mesh somewhere
  auto key = [shape = shape](const int& i, const int& j, const int& k,
                             const int& l) {
    return l + shape[3] * (k + shape[2] * (j + shape[1] * i));
  };

  int bin_key = key(i, j, k, l);

  if (bins.find(bin_key) == bins.end()) {
    bins[bin_key] = std::vector<BankedParticle*>();
  }

  bins[bin_key].push_back(&p);

  return true;
}

void ApproximateMeshCancelator::perform_cancellation(pcg32& /*rng*/) {
  // Go through all bins in the mesh
  for (auto& key_bin_pair : bins) {
    auto& bin = key_bin_pair.second;

    // Only do cancelation if we have more than one particle per bin
    if (bin.size() > 1) {
      // Only do cancellation if we have differing signs
      bool has_pos_w1 = false;
      bool has_neg_w1 = false;
      bool has_pos_w2 = false;
      bool has_neg_w2 = false;
      double sum_wgt = 0.;
      double sum_wgt2 = 0.;

      // Go through all particles in the bin
      for (const auto& p : bin) {
        if (p->wgt > 0.)
          has_pos_w1 = true;
        else if (p->wgt < 0.)
          has_neg_w1 = true;
        sum_wgt += p->wgt;

        if (p->wgt2 > 0.)
          has_pos_w2 = true;
        else if (p->wgt2 < 0.)
          has_neg_w2 = true;
        sum_wgt2 += p->wgt2;
      }

      // Get average weights
      double N = static_cast<double>(bin.size());
      double avg_wgt = sum_wgt / N;
      double avg_wgt2 = sum_wgt2 / N;

      // Go through all particles and change their weights
      for (auto& p : bin) {
        if (has_pos_w1 && has_neg_w1) p->wgt = avg_wgt;
        if (has_pos_w2 && has_neg_w2) p->wgt2 = avg_wgt2;
      }
    }

    bin.clear();
  }
}

std::vector<BankedParticle> ApproximateMeshCancelator::get_new_particles(
    pcg32& /*rng*/) {
  return {};
}

void ApproximateMeshCancelator::clear() { bins.clear(); }

std::shared_ptr<ApproximateMeshCancelator> make_approximate_mesh_cancelator(
    const YAML::Node& node) {
  // Get low
  if (!node["low"] || !node["low"].IsSequence() || !(node["low"].size() == 3)) {
    std::string mssg = "No valid low entry for approximate mesh cancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double xl = node["low"][0].as<double>();
  double yl = node["low"][1].as<double>();
  double zl = node["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!node["hi"] || !node["hi"].IsSequence() || !(node["hi"].size() == 3)) {
    std::string mssg = "No valid hi entry for approximate mesh cancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  double xh = node["hi"][0].as<double>();
  double yh = node["hi"][1].as<double>();
  double zh = node["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get shape
  if (!node["shape"] || !node["shape"].IsSequence() ||
      !(node["shape"].size() == 3)) {
    std::string mssg = "No valid shape entry for approximate mesh cancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  uint32_t Nx = node["shape"][0].as<uint32_t>();
  uint32_t Ny = node["shape"][1].as<uint32_t>();
  uint32_t Nz = node["shape"][2].as<uint32_t>();

  if (!node["energy-bounds"]) {
    return std::make_shared<ApproximateMeshCancelator>(r_low, r_hi, Nx, Ny, Nz);
  }

  if (!node["energy-bounds"].IsSequence()) {
    std::string mssg =
        "No valid energy-bounds entry for approximate mesh cancelator.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  std::vector<double> energy_bounds =
      node["energy-bounds"].as<std::vector<double>>();

  Output::instance()->write(" Using ApproximateMeshCancelator.\n");

  return std::make_shared<ApproximateMeshCancelator>(r_low, r_hi, Nx, Ny, Nz,
                                                     energy_bounds);
}
