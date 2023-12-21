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
#include <geometry/geometry.hpp>
#include <materials/material_helper.hpp>
#include <simulation/basic_exact_mg_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

#include <sobol/sobol.hpp>

#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

std::array<uint32_t, 3> BasicExactMGCancelator::KeyHash::shape;

BasicExactMGCancelator::BasicExactMGCancelator(Position low, Position hi,
                                               uint32_t Nx, uint32_t Ny,
                                               uint32_t Nz, BetaMode beta,
                                               bool sobol, uint32_t nsmp)
    : r_low(low),
      r_hi(hi),
      hash_fn(),
      dx((r_hi.x() - r_low.x()) / static_cast<double>(Nx)),
      dy((r_hi.y() - r_low.y()) / static_cast<double>(Ny)),
      dz((r_hi.z() - r_low.z()) / static_cast<double>(Nz)),
      beta_mode(beta),
      use_sobol(sobol),
      bins(),
      N_SAMPLES(nsmp) {
  // Make sure the low points are all lower than the high points
  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    fatal_error(
        "Low position is not lower than hi position in "
        "BasicExactMGCancelator.");
  }

  hash_fn.shape[0] = Nx;
  hash_fn.shape[1] = Ny;
  hash_fn.shape[2] = Nz;
}

bool BasicExactMGCancelator::add_particle(BankedParticle& p) {
  // Get bin indicies for spacial coordinates
  int i = static_cast<int>(std::floor((p.r.x() - r_low.x()) / dx));
  int j = static_cast<int>(std::floor((p.r.y() - r_low.y()) / dy));
  int k = static_cast<int>(std::floor((p.r.z() - r_low.z()) / dz));

  // If one of the indices is less than zero, don't keep the particle
  if (i < 0 || j < 0 || k < 0) {
    return false;
  }

  // if one of the indices is too large, don't keep the particle
  if (i >= static_cast<int>(hash_fn.shape[0]) ||
      j >= static_cast<int>(hash_fn.shape[1]) ||
      k >= static_cast<int>(hash_fn.shape[2])) {
    return false;
  }

  // Get the bin key
  Key key{i, j, k};

  // Get the material pointer
  Material* mat = get_material(p.r);

  // If the bin doesn't exist yet, initalize it
  if (bins.find(key) == bins.end()) {
    bins[key] = std::unordered_map<Material*, CancelBin>();
  }

  // Check if a bin exists for that material
  if (bins[key].find(mat) == bins[key].end()) {
    bins[key][mat] = CancelBin();
  }

  // Add particle and data
  bins[key][mat].particles.push_back(&p);
  bins[key][mat].W += p.wgt;
  bins[key][mat].W2 += p.wgt2;

  return true;
}

std::optional<Position> BasicExactMGCancelator::sample_position(
    const Key& key, Material* mat, pcg32& rng) const {
  // Get bin positions
  double Xl = r_low.x() + key.i * dx;
  double Yl = r_low.y() + key.j * dy;
  double Zl = r_low.z() + key.k * dz;

  uint32_t N_TRIES = 0;
  bool position_sampled = false;
  Position r_smp;
  while (N_TRIES < N_MAX_POS && !position_sampled) {
    double x = Xl + RNG::rand(rng) * dx;
    double y = Yl + RNG::rand(rng) * dy;
    double z = Zl + RNG::rand(rng) * dz;
    r_smp = Position(x, y, z);
    Material* mat_smp = get_material(r_smp);

    if (mat_smp == mat) {
      position_sampled = true;
    }

    N_TRIES++;
  }

  if (!position_sampled) {
    // We couldn't find a position with that material.
    return std::nullopt;
  }

  return std::make_optional(r_smp);
}

std::optional<Position> BasicExactMGCancelator::sample_position_sobol(
    const Key& key, Material* mat, unsigned long long& i) const {
  // Get bin positions
  double Xl = r_low.x() + key.i * dx;
  double Yl = r_low.y() + key.j * dy;
  double Zl = r_low.z() + key.k * dz;

  uint32_t N_TRIES = 0;
  bool position_sampled = false;
  Position r_smp;
  while (N_TRIES < N_MAX_POS && !position_sampled) {
    double x = Xl + sobol::sample(i, 0) * dx;
    double y = Yl + sobol::sample(i, 1) * dy;
    double z = Zl + sobol::sample(i, 2) * dz;
    r_smp = Position(x, y, z);
    Material* mat_smp = get_material(r_smp);

    if (mat_smp == mat) {
      position_sampled = true;
    }

    N_TRIES++;
    i++;
  }

  if (!position_sampled) {
    // We couldn't find a position with that material.
    return std::nullopt;
  }

  return std::make_optional(r_smp);
}

Material* BasicExactMGCancelator::get_material(const Position& r) const {
  Cell* cell = geometry::get_cell(r, {1., 0., 0.}).cell;

  if (!cell) {
    return nullptr;
  }

  Material* mat = cell->material();

  if (!mat) {
    std::stringstream mssg;
    mssg << "No material found at " << r
         << " in BasicExactMGCancelator::get_material.";
    fatal_error(mssg.str());
  }

  return mat;
}

double BasicExactMGCancelator::get_f(const Position& r,
                                     const Position& r_parent,
                                     double Esmp) const {
  double d = std::sqrt(std::pow(r.x() - r_parent.x(), 2.) +
                       std::pow(r.y() - r_parent.y(), 2.) +
                       std::pow(r.z() - r_parent.z(), 2.));

  // Can now calculate point fission reaction rate.
  // double f = (1. / (4. * PI * d * d)) * std::exp(-Esmp * d) * Ef;

  // In our [1], we used this commented version of f, with Ef and a
  // factor of 1/(4PI), but we have since learned that this is not necessary.
  // In simplified MG physics, where scattering is isotropic, and the fission
  // emisssion direction and energy are independent, and the fission energy is
  // independent of the incident energy, then we don't actually need to perfom
  // cancellation on the fission emission density, but only on the collision
  // density. As such, we can just use this modified version of f which is only
  // a function of the flight distance d, and the sampling xs Esmp. If the
  // fission specturm is a chi matrix however, or scattering is anisotropic,
  // this will not work !
  // We will leave this for now with the compiler warning due to the unused Ef,
  // as this must be changed to accomodate the new MG physics which allows for
  // anisotropic scattering, and chi matrix.
  double f = (1. / (d * d)) * std::exp(-Esmp * d);

  return f;
}

double BasicExactMGCancelator::get_min_f(const Key& key,
                                         const Position& r_parent,
                                         double Esmp) const {
  // Get bin positions
  double Xl = r_low.x() + key.i * dx;
  double Xh = Xl + dx;
  double Yl = r_low.y() + key.j * dy;
  double Yh = Yl + dy;
  double Zl = r_low.z() + key.k * dz;
  double Zh = Zl + dz;

  // Create 8 positions for cube corners
  Position c1(Xh, Yh, Zh);
  Position c2(Xh, Yh, Zl);
  Position c3(Xh, Yl, Zh);
  Position c4(Xh, Yl, Zl);
  Position c5(Xl, Yh, Zh);
  Position c6(Xl, Yh, Zl);
  Position c7(Xl, Yl, Zh);
  Position c8(Xl, Yl, Zl);

  double f1 = get_f(c1, r_parent, Esmp);
  double f2 = get_f(c2, r_parent, Esmp);
  double f3 = get_f(c3, r_parent, Esmp);
  double f4 = get_f(c4, r_parent, Esmp);
  double f5 = get_f(c5, r_parent, Esmp);
  double f6 = get_f(c6, r_parent, Esmp);
  double f7 = get_f(c7, r_parent, Esmp);
  double f8 = get_f(c8, r_parent, Esmp);

  double beta = std::min(std::min(std::min(f1, f2), std::min(f3, f4)),
                         std::min(std::min(f5, f6), std::min(f7, f8)));

  return beta;
}

void BasicExactMGCancelator::get_averages(const Key& key, Material* mat,
                                          CancelBin& bin, pcg32& rng) {
  // Make sure vectors are allocated
  bin.averages.resize(bin.particles.size());

  // Get all positions
  std::vector<Position> r_smps;
  r_smps.reserve(N_SAMPLES);
  for (std::size_t j = 0; j < N_SAMPLES; j++) {
    std::optional<Position> r_smp = sample_position(key, mat, rng);
    if (r_smp) {
      r_smps.push_back(r_smp.value());
    } else {
      // We couldn't sample a point, so we just wont cancel
      // this bin.
      bin.can_cancel = false;
      return;
    }
  }

  // Go through all particles
  for (std::size_t i = 0; i < bin.particles.size(); i++) {
    Position r_parent = bin.particles[i]->parents_previous_position;
    double E = bin.particles[i]->parents_previous_energy;
    double Esmp = bin.particles[i]->Esmp_parent;
    MaterialHelper mat_helper(mat, E);

    // Compute the average value for f and 1/f
    double sum_f = 0.;
    double sum_f_inv = 0.;
    for (const auto& r_smp : r_smps) {
      double f = get_f(r_smp, r_parent, Esmp);
      sum_f += f;
      sum_f_inv += 1. / f;
    }
    bin.averages[i].f = sum_f / static_cast<double>(N_SAMPLES);
    bin.averages[i].f_inv = sum_f_inv / static_cast<double>(N_SAMPLES);
  }

  if (beta_mode == BetaMode::OptAverageGain) {
    auto C = [](double f, double f_inv) { return 1. / (2. * f * f_inv - 1.); };

    double sum_c = 0.;
    for (std::size_t i = 0; i < bin.particles.size(); i++)
      sum_c += C(bin.averages[i].f, bin.averages[i].f_inv);

    double sum_c_wgt = 0.;
    double sum_c_wgt2 = 0.;
    for (std::size_t i = 0; i < bin.particles.size(); i++) {
      sum_c_wgt +=
          C(bin.averages[i].f, bin.averages[i].f_inv) * bin.particles[i]->wgt;
      sum_c_wgt2 +=
          C(bin.averages[i].f, bin.averages[i].f_inv) * bin.particles[i]->wgt2;
    }

    bin.sum_c = sum_c;
    bin.sum_c_wgt = sum_c_wgt;
    bin.sum_c_wgt2 = sum_c_wgt2;
  }
}

void BasicExactMGCancelator::get_averages_sobol(const Key& key, Material* mat,
                                                CancelBin& bin) {
  // Make sure vectors are allocated
  bin.averages.resize(bin.particles.size());

  // Get all positions
  std::vector<Position> r_smps;
  r_smps.reserve(N_SAMPLES);
  unsigned long long sobol_index = 0;
  for (std::size_t j = 0; j < N_SAMPLES; j++) {
    std::optional<Position> r_smp =
        sample_position_sobol(key, mat, sobol_index);
    if (r_smp) {
      r_smps.push_back(r_smp.value());
    } else {
      // We couldn't sample a point, so we just wont cancel
      // this bin.
      bin.can_cancel = false;
      return;
    }
  }

  // Go through all particles
  for (std::size_t i = 0; i < bin.particles.size(); i++) {
    Position r_parent = bin.particles[i]->parents_previous_position;
    double E = bin.particles[i]->parents_previous_energy;
    double Esmp = bin.particles[i]->Esmp_parent;
    MaterialHelper mat_helper(mat, E);

    // Compute the average value for f and 1/f
    double sum_f = 0.;
    double sum_f_inv = 0.;
    for (const auto& r_smp : r_smps) {
      double f = get_f(r_smp, r_parent, Esmp);
      sum_f += f;
      sum_f_inv += 1. / f;
    }
    bin.averages[i].f = sum_f / static_cast<double>(N_SAMPLES);
    bin.averages[i].f_inv = sum_f_inv / static_cast<double>(N_SAMPLES);
  }

  if (beta_mode == BetaMode::OptAverageGain) {
    auto C = [](double f, double f_inv) { return 1. / (2. * f * f_inv - 1.); };

    double sum_c = 0.;
    for (std::size_t i = 0; i < bin.particles.size(); i++)
      sum_c += C(bin.averages[i].f, bin.averages[i].f_inv);

    double sum_c_wgt = 0.;
    double sum_c_wgt2 = 0.;
    for (std::size_t i = 0; i < bin.particles.size(); i++) {
      sum_c_wgt +=
          C(bin.averages[i].f, bin.averages[i].f_inv) * bin.particles[i]->wgt;
      sum_c_wgt2 +=
          C(bin.averages[i].f, bin.averages[i].f_inv) * bin.particles[i]->wgt2;
    }

    bin.sum_c = sum_c;
    bin.sum_c_wgt = sum_c_wgt;
    bin.sum_c_wgt2 = sum_c_wgt2;
  }
}

double BasicExactMGCancelator::get_beta(const Key& key, const CancelBin& bin,
                                        std::size_t i, const Position& r_parent,
                                        double Esmp, double wgt,
                                        bool first_wgt) const {
  if (bin.can_cancel == false) return 0.;

  switch (beta_mode) {
    case BetaMode::Zero:
      return 0.;
      break;

    case BetaMode::Minimum:
      return get_min_f(key, r_parent, Esmp);
      break;

    case BetaMode::OptAverageF: {
      double f = bin.averages[i].f;
      double N = static_cast<double>(bin.particles.size());
      double W = first_wgt ? bin.W : bin.W2;
      return f * (1. - ((W) / ((N + 1.) * wgt)));
      break;
    }

    case BetaMode::OptAverageGain: {
      double sum_c_wgt = first_wgt ? bin.sum_c_wgt : bin.sum_c_wgt2;
      double S = sum_c_wgt / (1. + bin.sum_c);
      double f = bin.averages[i].f;
      double f_inv = bin.averages[i].f_inv;
      return f * (1. / (2. * f * f_inv - 1.)) * (1. - (S / wgt));
      break;
    }
  }

  // Never gets here !
  return 0.;
}

void BasicExactMGCancelator::cancel_bin(const Key& key, Material* mat_ptr,
                                        CancelBin& bin, bool first_wgt) {
  // Go through all particles in bin
  for (std::size_t i = 0; i < bin.particles.size(); i++) {
    Position r_parent = bin.particles[i]->parents_previous_position;
    double E = bin.particles[i]->parents_previous_energy;
    double Esmp = bin.particles[i]->Esmp_parent;
    MaterialHelper mat(mat_ptr, E);

    double wgt = first_wgt ? bin.particles[i]->wgt : bin.particles[i]->wgt2;

    // If the weight component is zero, we can't cancel that part.
    // This is required as one component can be zero in noise simulations.
    if (wgt == 0.) return;

    const double B = get_beta(key, bin, i, r_parent, Esmp, wgt, first_wgt);
    const double f = get_f(bin.particles[i]->r, r_parent, Esmp);

    const double P_p = (f - B) / f;
    const double P_u = B / f;

    // If either P_p or P_u is Inf or NaN, we can't do cancellation.
    // This usually happens when the weight in question is zero, which
    // occurs durring noise transport. This should be checked for just
    // above, but I am gonna leave this here anyway for safety.
    if (std::isinf(P_u) || std::isinf(P_p) || std::isnan(P_u) ||
        std::isnan(P_p))
      return;

    if (first_wgt) {
      bin.uniform_wgt += bin.particles[i]->wgt * P_u;
      bin.particles[i]->wgt *= P_p;
    } else {
      bin.uniform_wgt2 += bin.particles[i]->wgt2 * P_u;
      bin.particles[i]->wgt2 *= P_p;
    }
  }
}

void BasicExactMGCancelator::perform_cancellation(pcg32& rng) {
  if (beta_mode == BetaMode::Zero) return;

  // If we have no bins (meaning no particles), then
  // we can't do any cancellation.
  if (bins.size() == 0) return;

  // Get vector of keys to do cancellation in parallel
  std::vector<std::pair<Key, Material*>> keys;
  keys.reserve(bins.size());
  for (const auto& key_matbin_pair : bins) {
    const auto& key = key_matbin_pair.first;
    const auto& matbin = key_matbin_pair.second;
    for (const auto& mat_bin_pair : matbin)
      keys.push_back({key, mat_bin_pair.first});
  }

  // Get seed offsets to try and make parallel cancellation
  // deterministic (i.e. independent of the number of threads)
  uint64_t seed_advance = 0;
  if (beta_mode != BetaMode::Minimum) {
    uint64_t max_rn_per_part =
        N_SAMPLES * N_MAX_POS * (beta_mode == BetaMode::OptAverageGain ? 2 : 1);
    bins[keys[0].first][keys[0].second].rng_seed_advance = 0;
    seed_advance =
        bins[keys[0].first][keys[0].second].particles.size() * max_rn_per_part;
    for (std::size_t i = 1; i < keys.size(); i++) {
      bins[keys[i].first][keys[i].second].rng_seed_advance = seed_advance;
      seed_advance += bins[keys[i].first][keys[i].second].particles.size() *
                      max_rn_per_part;
    }
  }

  // Go through all bins
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < keys.size(); i++) {
    Key key = keys[i].first;
    Material* mat = keys[i].second;
    CancelBin& bin = bins[key][mat];

    pcg32 rng_local = rng;
    rng_local.advance(bin.rng_seed_advance);

    // Only atempt cancelation if we have two or more particles
    if (bin.particles.size() > 1) {
      // Go through all particles and see if we have mixed signes
      bool has_pos_w1 = false;
      bool has_neg_w1 = false;
      bool has_pos_w2 = false;
      bool has_neg_w2 = false;
      for (const auto& p : bin.particles) {
        if (p->wgt > 0.)
          has_pos_w1 = true;
        else if (p->wgt < 0.)
          has_neg_w1 = true;

        if (p->wgt2 > 0.)
          has_pos_w2 = true;
        else if (p->wgt2 < 0.)
          has_neg_w2 = true;

        if (has_pos_w1 && has_neg_w1 && has_pos_w2 && has_neg_w2) break;
      }

      // If we are doing cancellation and the beta_mode requires it, get the
      // average value of the fission density and 1 / fission density for each
      // particle in the bin.
      if (((has_pos_w1 && has_neg_w1) || (has_pos_w2 && has_neg_w2)) &&
          (beta_mode == BetaMode::OptAverageF ||
           beta_mode == BetaMode::OptAverageGain)) {
        if (use_sobol) {
          get_averages_sobol(key, mat, bin);
        } else {
          get_averages(key, mat, bin, rng_local);
        }
      }

      // If determined necessary, carry out the cancellations for each weight
      if (has_pos_w1 && has_neg_w1) cancel_bin(key, mat, bin, true);
      if (has_pos_w2 && has_neg_w2) cancel_bin(key, mat, bin, false);

      // Cancellation has now occured. We should clear the bin of particles
      bin.particles.clear();
      bin.averages.clear();
      bin.sum_c = 0.;
      bin.sum_c_wgt = 0.;
      bin.sum_c_wgt2 = 0.;
    }
  }

  if (beta_mode != BetaMode::Minimum) {
    rng.advance(seed_advance);
  }
}

std::vector<BankedParticle> BasicExactMGCancelator::get_new_particles(
    pcg32& rng) {
  if (beta_mode == BetaMode::Zero) return {};

  std::vector<BankedParticle> uniform_particles;

  for (auto& key_bin_pair : bins) {
    const auto& key = key_bin_pair.first;
    auto& material_bins = key_bin_pair.second;

    for (auto& mat_bin_pair : material_bins) {
      Material* mat = mat_bin_pair.first;
      CancelBin& bin = mat_bin_pair.second;

      // Get bin positions
      auto mat_ptr = mat->shared_from_this();

      // Determine number of new particles to add
      uint32_t N = static_cast<uint32_t>(std::ceil(
          std::max(std::abs(bin.uniform_wgt), std::abs(bin.uniform_wgt2))));

      if (N > 0) {
        double w = bin.uniform_wgt / N;
        double w2 = bin.uniform_wgt2 / N;

        // Get the single nuclide from the material
        std::shared_ptr<Nuclide> nuclide =
            mat_ptr->components().front().nuclide;

        for (size_t i = 0; i < N; i++) {
          // Sample position
          std::optional<Position> r_smp = sample_position(key, mat, rng);

          if (!r_smp) {
            // For some reason we couldn't sample a position, probably because
            // there is so little of the material in the region. If this
            // is the case, we shouldn't have gotten this far, but hey, here
            // were are ! We are just gonna call this a fatal error for now,
            // and see if it ever pops up.
            fatal_error("Couldn't sample position for uniform particle.");
          }

          // Sample particle energy and direction
          FissionInfo finfo =
              nuclide->sample_fission(0., {0., 0., 1.}, 0, 0., rng);
          BankedParticle uniform_particle;
          uniform_particle.r = r_smp.value();
          uniform_particle.u = finfo.direction;
          uniform_particle.E = finfo.energy;
          uniform_particle.wgt = w;
          uniform_particle.wgt2 = w2;
          uniform_particle.parent_history_id = 0;
          uniform_particle.parent_daughter_id = 0;

          // Save sampled particle
          uniform_particles.push_back(uniform_particle);
        }
      }

      bin.uniform_wgt = 0.;
      bin.uniform_wgt2 = 0.;
    }
  }

  return uniform_particles;
}

void BasicExactMGCancelator::clear() { bins.clear(); }

std::shared_ptr<BasicExactMGCancelator> make_basic_exact_mg_cancelator(
    const YAML::Node& node) {
  // Get low
  if (!node["low"] || !node["low"].IsSequence() || !(node["low"].size() == 3)) {
    fatal_error("No valid low entry for basic exact MG cancelator.");
  }

  double xl = node["low"][0].as<double>();
  double yl = node["low"][1].as<double>();
  double zl = node["low"][2].as<double>();

  Position r_low(xl, yl, zl);

  // Get hi
  if (!node["hi"] || !node["hi"].IsSequence() || !(node["hi"].size() == 3)) {
    fatal_error("No valid hi entry for basic exact MG cancelator.");
  }

  double xh = node["hi"][0].as<double>();
  double yh = node["hi"][1].as<double>();
  double zh = node["hi"][2].as<double>();

  Position r_hi(xh, yh, zh);

  // Get shape
  if (!node["shape"] || !node["shape"].IsSequence() ||
      !(node["shape"].size() == 3)) {
    fatal_error("No valid shape entry for basic exact MG cancelator.");
  }

  uint32_t Nx = node["shape"][0].as<uint32_t>();
  uint32_t Ny = node["shape"][1].as<uint32_t>();
  uint32_t Nz = node["shape"][2].as<uint32_t>();

  if (!node["beta"] || !node["beta"].IsScalar()) {
    fatal_error("No valid beta entry for basic exact MG cancelator.");
  }
  std::string beta_str = node["beta"].as<std::string>();

  BasicExactMGCancelator::BetaMode beta =
      BasicExactMGCancelator::BetaMode::Zero;

  if (beta_str == "zero") {
    beta = BasicExactMGCancelator::BetaMode::Zero;
  } else if (beta_str == "minimum") {
    beta = BasicExactMGCancelator::BetaMode::Minimum;
  } else if (beta_str == "average-f") {
    beta = BasicExactMGCancelator::BetaMode::OptAverageF;
  } else if (beta_str == "average-g") {
    beta = BasicExactMGCancelator::BetaMode::OptAverageGain;
  } else {
    fatal_error("Unkown beta entry \"" + beta_str +
                "\" for basic exact MG cancelator.");
  }

  bool use_sobol = true;
  if (node["sobol"] && node["sobol"].IsScalar()) {
    use_sobol = node["sobol"].as<bool>();
  } else if (node["sobol"]) {
    fatal_error("Invalid sobol entry for cancelator.");
  }

  uint32_t n_samples = 10;
  if (node["n-samples"] && node["n-samples"].IsScalar()) {
    n_samples = node["n-samples"].as<uint32_t>();
  } else if (node["n-samples"]) {
    fatal_error("Invalid n-samples entry for cancelator.");
  }

  if (n_samples == 0) {
    fatal_error("n-samples must be greater than zero.");
  }

  std::stringstream otpt;
  otpt << " Using BasicExactMGCancelator with BetaMode " << beta_str;
  if (use_sobol) otpt << " with sobol";
  if (beta != BasicExactMGCancelator::BetaMode::Zero &&
      beta != BasicExactMGCancelator::BetaMode::Minimum) {
    otpt << " using " << n_samples << " points";
  }
  otpt << ".\n";
  Output::instance().write(otpt.str());

  return std::make_shared<BasicExactMGCancelator>(r_low, r_hi, Nx, Ny, Nz, beta,
                                                  use_sobol, n_samples);
}

//==============================================================================
// References
//
// [1] H. Belanger, D. Mancusi, and A. Zoia, “Exact weight cancellation in Monte
//     Carlo eigenvalue transport problems,” Phys. Rev. E, vol. 104, no. 1,
//     p. 015306, 2021, doi: 10.1103/physreve.104.015306.
