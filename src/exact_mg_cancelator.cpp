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
#include <simulation/exact_mg_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>

#include <sobol/sobol.hpp>

#include <cmath>
#include <set>

Position ExactMGCancelator::Key::r_low, ExactMGCancelator::Key::r_hi;
std::array<std::size_t, 4> ExactMGCancelator::Key::shape;
std::array<double, 3> ExactMGCancelator::Key::pitch;
std::vector<std::vector<std::size_t>> ExactMGCancelator::Key::group_bins;

ExactMGCancelator::ExactMGCancelator(
    const Position& r_low, const Position& r_hi,
    const std::array<std::size_t, 4>& shape,
    const std::vector<std::vector<std::size_t>>& group_bins, bool chi_matrix,
    bool use_virtual_collisions, uint32_t n_samples)
    : bins(),
      CHI_MATRIX(chi_matrix),
      N_SAMPLES(n_samples),
      USE_VIRTUAL_COLLISIONS(use_virtual_collisions) {
  // Make sure N_SAMPLES != 0
  if (N_SAMPLES == 0) {
    fatal_error("N_SAMPLES == 0");
  }

  // Make sure r_low is low
  if (r_low.x() >= r_hi.x()) {
    fatal_error("r_low.x() >= r_hi.x()");
  }
  if (r_low.y() >= r_hi.y()) {
    fatal_error("r_low.y() >= r_hi.y()");
  }
  if (r_low.z() >= r_hi.z()) {
    fatal_error("r_low.z() >= r_hi.z()");
  }

  // Make sure shape[3] == group_bins.size()
  if (shape[3] != group_bins.size()) {
    fatal_error("shape[3] != group_bins.size()");
  }

  if (CHI_MATRIX && shape[3] == 0) {
    fatal_error(
        "Chi matrix is used, but no group_bins provided.\nImpossible to have "
        "exact cancellation.");
  }

  std::set<std::size_t> bined_groups;
  for (const auto& bin : group_bins) {
    for (const auto& grp : bin) {
      if (grp >= settings::ngroups) {
        fatal_error("grp >= settings::ngroups");
      }

      if (bined_groups.count(grp) == 1) {
        std::stringstream mssg;
        mssg << "Group " << grp << " is included more than once in group_bins.";
        fatal_error(mssg.str());
      } else {
        bined_groups.insert(grp);
      }
    }
  }

  // Set the key values
  Key::r_low = r_low;
  Key::r_hi = r_hi;
  Key::shape[0] = shape[0];
  Key::shape[1] = shape[1];
  Key::shape[2] = shape[2];
  Key::shape[3] = shape[3];
  Key::pitch[0] = (r_hi.x() - r_low.x()) / static_cast<double>(shape[0]);
  Key::pitch[1] = (r_hi.y() - r_low.y()) / static_cast<double>(shape[1]);
  Key::pitch[2] = (r_hi.z() - r_low.z()) / static_cast<double>(shape[2]);
  Key::group_bins = group_bins;

  // If we have CHI_MATRIX == false, then we do not need to
  // provide any group-bins to the cancelator to have exact
  // cancellation. However, if this is the case, the hash
  // function will always return 0 if shape[3] == 0, and
  // we will get too many collisions. As such, we must set
  // it equal to 1.
  if (Key::shape[3] == 0) Key::shape[3] = 1;
}

std::optional<ExactMGCancelator::Key> ExactMGCancelator::get_key(
    const Position& r, std::size_t g) {
  // If we are outside the spatial mesh, return nullopt.
  if (r.x() < Key::r_low.x() || r.x() > Key::r_hi.x() ||
      r.y() < Key::r_low.y() || r.y() > Key::r_hi.y() ||
      r.z() < Key::r_low.z() || r.z() > Key::r_hi.z()) {
    return std::nullopt;
  }

  // Get the spatial indices.
  std::size_t i = static_cast<std::size_t>(
      std::floor((r.x() - Key::r_low.x()) / Key::pitch[0]));
  std::size_t j = static_cast<std::size_t>(
      std::floor((r.y() - Key::r_low.y()) / Key::pitch[1]));
  std::size_t k = static_cast<std::size_t>(
      std::floor((r.z() - Key::r_low.z()) / Key::pitch[2]));

  // Now we need the energy index.
  bool e_determined = false;
  std::size_t e = 0;

  // Go through all energy bins
  for (e = 0; e < Key::group_bins.size(); e++) {
    for (std::size_t i = 0; i < Key::group_bins[e].size(); i++) {
      if (g == Key::group_bins[e][i]) {
        e_determined = true;
        break;
      }
    }
    if (e_determined) break;
  }

  // Make sure we have a valid energy bin index
  if (e_determined == false && CHI_MATRIX) return std::nullopt;

  // Return the key
  return std::make_optional<Key>(Key(i, j, k, e));
}

Material* ExactMGCancelator::get_material(const Position& r) const {
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

bool ExactMGCancelator::add_particle(BankedParticle& p) {
  // Get the key for p
  std::size_t g_fiss = settings::group(p.E);
  auto optional_key = get_key(p.r, g_fiss);

  // If we don't get a key, then the particle isn't inside
  // any cancellation region.
  if (!optional_key) return false;

  auto key = *optional_key;

  if (USE_VIRTUAL_COLLISIONS == false && p.parents_previous_was_virtual) {
    // If the parent's previous collision was virtual, we cannot reliably
    // perform cancellation in an accurate manner. This is due to the
    // fact that we can't scatter forward, and not change energy.
    // As such, we do not perform cancellation on particles which were
    // produced imediately after a virtual collision. We return true
    // however, as the particle was placed in a valid bin.
    return true;
  }

  // If the bin doesn't exist yet, initalize it
  if (bins.find(key) == bins.end()) {
    bins[key] = std::unordered_map<Material*, CancelBin>();
  }

  // Get the pointer to the material where the particle is
  Material* mat = get_material(p.r);

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

double ExactMGCancelator::get_f(const Position& r1, const Direction& u1,
                                std::size_t g1, std::size_t g3,
                                const Position& r4, std::size_t g4, double Esmp,
                                MGNuclide* nuclide) const {
  // Get distance between r1 an r4
  const double d = (r4 - r1).norm();

  // Direction from r1 to r4
  const Direction u(r4.x() - r1.x(), r4.y() - r1.y(), r4.z() - r1.z());

  // Get scattering cosine
  const double mu = u.dot(u1);

  // Get scattering pdf
  const double pdf_mu = nuclide->angles()[g1][g3].pdf(mu);

  // Get the pdf for the chi portion
  const double pdf_chi = CHI_MATRIX ? nuclide->chi()[g3][g4] : 1.;

  return (pdf_mu * pdf_chi / (d * d)) * std::exp(-Esmp * d);
}

double ExactMGCancelator::get_beta(const CancelBin& bin, std::size_t i,
                                   bool wgt_1) const {
  if (bin.can_cancel == false) return 0.;

  const double wgt = wgt_1 ? bin.particles[i]->wgt : bin.particles[i]->wgt2;
  const double sum_c_wgt = wgt_1 ? bin.sum_c_wgt : bin.sum_c_wgt2;
  const double S = sum_c_wgt / (1. + bin.sum_c);
  const double f = bin.averages[i].f;
  const double f_inv = bin.averages[i].f_inv;
  return f * (1. / (2. * f * f_inv - 1.)) * (1. - (S / wgt));
}

std::optional<std::pair<Position, std::size_t>> ExactMGCancelator::sample_point(
    const Key& key, Material* mat, unsigned long long& i) const {
  // Get lower limits of bin
  double Xl = Key::r_low.x() + static_cast<double>(key.i) * Key::pitch[0];
  double Yl = Key::r_low.y() + static_cast<double>(key.j) * Key::pitch[1];
  double Zl = Key::r_low.z() + static_cast<double>(key.k) * Key::pitch[2];

  uint32_t N_TRIES = 0;
  bool position_sampled = false;
  Position r_smp;
  while (N_TRIES < N_MAX_POS && !position_sampled) {
    double x = Xl + sobol::sample(i, 0) * Key::pitch[0];
    double y = Yl + sobol::sample(i, 1) * Key::pitch[1];
    double z = Zl + sobol::sample(i, 2) * Key::pitch[2];
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

  // Fission group
  std::size_t g = 0;

  // Now that we have a position, we need to sample
  // an energy group if we have chi matrices in the problem.
  if (CHI_MATRIX) {
    // Get the 4th sobol coordinate for i-1, as we used i-1
    // to sample the position (loop above does i++ at end).
    const double xi = sobol::sample(i - 1, 3);
    const std::size_t g_index = static_cast<std::size_t>(
        std::floor(xi * static_cast<double>(Key::group_bins[key.e].size())));
    g = Key::group_bins[key.e][g_index];
  }

  return std::make_optional<std::pair<Position, std::size_t>>({r_smp, g});
}

void ExactMGCancelator::compute_averages(const Key& key, Material* mat,
                                         MGNuclide* nuclide, CancelBin& bin) {
  // Make sure vectors are allocated
  bin.averages.resize(bin.particles.size());

  // Get all position and energies
  std::vector<std::pair<Position, std::size_t>> r_E_smps;
  r_E_smps.reserve(N_SAMPLES);
  unsigned long long sobol_index = 0;
  for (std::size_t j = 0; j < N_SAMPLES; j++) {
    std::optional<std::pair<Position, std::size_t>> r_E_smp =
        sample_point(key, mat, sobol_index);
    if (r_E_smp) {
      r_E_smps.push_back(r_E_smp.value());
    } else {
      // We couldn't sample a point, so we just wont cancel
      // this bin.
      bin.can_cancel = false;
      return;
    }
  }

  // Go through all particles
  for (std::size_t i = 0; i < bin.particles.size(); i++) {
    const Position& r1 = bin.particles[i]->parents_previous_position;
    const Direction& u1 = bin.particles[i]->parents_previous_direction;
    const double E1 = bin.particles[i]->parents_previous_previous_energy;
    const std::size_t g1 = settings::group(E1);
    const double E3 = bin.particles[i]->parents_previous_energy;
    const std::size_t g3 = settings::group(E3);
    const double Esmp = bin.particles[i]->Esmp_parent;

    // Compute the average value for f and 1/f
    double sum_f = 0.;
    double sum_f_inv = 0.;
    for (const auto& r_E_smp : r_E_smps) {
      const Position& r4 = r_E_smp.first;
      const std::size_t g4 = r_E_smp.second;
      double f = get_f(r1, u1, g1, g3, r4, g4, Esmp, nuclide);
      if (f == 0.) {
        bin.can_cancel = false;
        return;
      }
      sum_f += f;
      sum_f_inv += 1. / f;
    }
    bin.averages[i].f = sum_f / static_cast<double>(N_SAMPLES);
    bin.averages[i].f_inv = sum_f_inv / static_cast<double>(N_SAMPLES);
  }

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

void ExactMGCancelator::cancel_bin(CancelBin& bin, MGNuclide* nuclide,
                                   bool first_wgt) {
  // Go through all particles in bin
  for (std::size_t i = 0; i < bin.particles.size(); i++) {
    const Position& r1 = bin.particles[i]->parents_previous_position;
    const Direction& u1 = bin.particles[i]->parents_previous_direction;
    const double E1 = bin.particles[i]->parents_previous_previous_energy;
    const std::size_t g1 = settings::group(E1);
    const double E3 = bin.particles[i]->parents_previous_energy;
    const std::size_t g3 = settings::group(E3);
    const std::size_t g4 = settings::group(bin.particles[i]->E);
    const Position& r4 = bin.particles[i]->r;
    const double Esmp = bin.particles[i]->Esmp_parent;

    const double B = get_beta(bin, i, first_wgt);
    const double f = get_f(r1, u1, g1, g3, r4, g4, Esmp, nuclide);

    const double P_p = (f - B) / f;
    const double P_u = B / f;

    // If either P_p or P_u is Inf or NaN, we can't do cancellation.
    // This usually happens when the weight in question is zero, which
    // occurs durring noise transport.
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

void ExactMGCancelator::perform_cancellation(pcg32&) {
  // If we have no bins (meaning no particles), then
  // we can't do any cancellation.
  if (bins.size() == 0) return;

  // Get vector of keys to do cancellation in parallel
  std::vector<std::pair<Key, Material*>> key_mat_pairs;
  key_mat_pairs.reserve(bins.size());
  for (const auto& key_matbin_pair : bins) {
    const auto& key = key_matbin_pair.first;
    const auto& matbin = key_matbin_pair.second;
    for (const auto& mat_bin_pair : matbin)
      key_mat_pairs.emplace_back(key, mat_bin_pair.first);
  }

  // Go through all bins
#ifdef ABEILLE_USE_OMP
#pragma omp parallel for
#endif
  for (std::size_t i = 0; i < key_mat_pairs.size(); i++) {
    Key key = key_mat_pairs[i].first;
    Material* mat = key_mat_pairs[i].second;
    CancelBin& bin = bins[key][mat];

    // For cancellation to be eact in the most general MG case,
    // need access to all scattering PDFs, and to the Chi matrix,
    // if present. Since ExactMGCancelator should only be used in
    // MG mode, this static cast should be safe, as no CENuclide
    // instances should have been created. This check is done in
    // make_cancelator in cancellator.cpp.
    MGNuclide* nuclide =
        static_cast<MGNuclide*>(mat->components()[0].nuclide.get());

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
      if ((has_pos_w1 && has_neg_w1) || (has_pos_w2 && has_neg_w2)) {
        compute_averages(key, mat, nuclide, bin);
      }

      // If determined necessary, carry out the cancellations for each weight
      if (has_pos_w1 && has_neg_w1) cancel_bin(bin, nuclide, true);
      if (has_pos_w2 && has_neg_w2) cancel_bin(bin, nuclide, false);

      // Cancellation has now occured. We should clear the bin of particles
      bin.particles.clear();
      bin.averages.clear();
      bin.sum_c = 0.;
      bin.sum_c_wgt = 0.;
      bin.sum_c_wgt2 = 0.;
    }
  }
}

std::optional<Position> ExactMGCancelator::sample_position(const Key& key,
                                                           Material* mat,
                                                           pcg32& rng) const {
  // Get bin positions
  double Xl = Key::r_low.x() + static_cast<double>(key.i) * Key::pitch[0];
  double Yl = Key::r_low.y() + static_cast<double>(key.j) * Key::pitch[1];
  double Zl = Key::r_low.z() + static_cast<double>(key.k) * Key::pitch[2];

  uint32_t N_TRIES = 0;
  bool position_sampled = false;
  Position r_smp;
  while (N_TRIES < N_MAX_POS && !position_sampled) {
    double x = Xl + RNG::rand(rng) * Key::pitch[0];
    double y = Yl + RNG::rand(rng) * Key::pitch[1];
    double z = Zl + RNG::rand(rng) * Key::pitch[2];
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

std::vector<BankedParticle> ExactMGCancelator::get_new_particles(pcg32& rng) {
  std::vector<BankedParticle> uniform_particles;

  for (auto& key_bin_pair : bins) {
    const auto& key = key_bin_pair.first;
    auto& material_bins = key_bin_pair.second;

    for (auto& mat_bin_pair : material_bins) {
      Material* mat = mat_bin_pair.first;
      CancelBin& bin = mat_bin_pair.second;

      // Determine number of new particles to add
      uint32_t N = static_cast<uint32_t>(std::ceil(
          std::max(std::abs(bin.uniform_wgt), std::abs(bin.uniform_wgt2))));

      if (N > 0) {
        double w = bin.uniform_wgt / N;
        double w2 = bin.uniform_wgt2 / N;

        // Get the single nuclide from the material, and convert to MGNuclide.
        // This SHOULD be safe, as in theory, we wont construct an
        // ExactMGCancelator instance in CE mode.
        MGNuclide* nuclide =
            static_cast<MGNuclide*>(mat->components()[0].nuclide.get());

        // Make all N uniform particles for this bin
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

          // Now we sample an energy group
          std::size_t e_index = 0;
          if (CHI_MATRIX) {
            // We need to select a random energy group from our bin.
            const double xi_E = RNG::rand(rng);
            std::size_t g_index = static_cast<std::size_t>(std::floor(
                xi_E * static_cast<double>(Key::group_bins[key.e].size())));
            e_index = Key::group_bins[key.e][g_index];
          } else {
            // Fission spectrum is independent of incident energy.
            // We can just sample an energy from the first row
            // of the chi matrix.
            e_index =
                static_cast<std::size_t>(RNG::discrete(rng, nuclide->chi()[0]));
          }
          double E_smp = 0.5 * (settings::energy_bounds[e_index] +
                                settings::energy_bounds[e_index + 1]);

          // Sample Direction
          Direction u_smp(2. * RNG::rand(rng) - 1., 2. * PI * RNG::rand(rng));

          // Construct uniform_particle
          BankedParticle uniform_particle;
          uniform_particle.r = r_smp.value();
          uniform_particle.u = u_smp;
          uniform_particle.E = E_smp;
          uniform_particle.wgt = w;
          uniform_particle.wgt2 = w2;

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

void ExactMGCancelator::clear() { bins.clear(); }

std::shared_ptr<ExactMGCancelator> make_exact_mg_cancelator(
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

  if (r_low.x() >= r_hi.x() || r_low.y() >= r_hi.y() || r_low.z() >= r_hi.z()) {
    fatal_error("In ExactMGCancelator, low is greater than or equal to hi.");
  }

  // Get shape
  if (!node["shape"] || !node["shape"].IsSequence() ||
      !(node["shape"].size() == 3)) {
    fatal_error("No valid shape entry for basic exact MG cancelator.");
  }

  uint32_t Nx = node["shape"][0].as<uint32_t>();
  uint32_t Ny = node["shape"][1].as<uint32_t>();
  uint32_t Nz = node["shape"][2].as<uint32_t>();
  uint32_t Ne = 0;

  std::vector<std::vector<std::size_t>> group_bins;
  if (node["group-bins"] && node["group-bins"].IsSequence()) {
    group_bins = node["group-bins"].as<std::vector<std::vector<std::size_t>>>();
  } else if (settings::chi_matrix) {
    warning(
        "No group-bins are provided to ExactMGCancelator, but fission spectra "
        "are provided with a matrix.\nCancellation cannot be exact.");
  }
  Ne = static_cast<uint32_t>(group_bins.size());

  std::array<std::size_t, 4> shape{Nx, Ny, Nz, Ne};

  std::set<std::size_t> bined_groups;
  for (const auto& bin : group_bins) {
    for (const auto& grp : bin) {
      if (grp >= settings::ngroups) {
        std::stringstream mssg;
        mssg << "Invalid group " << grp << " in group-bins.";
        fatal_error(mssg.str());
      }

      if (bined_groups.count(grp) == 1) {
        std::stringstream mssg;
        mssg << "Group " << grp << " is included more than once in group-bins.";
        fatal_error(mssg.str());
      } else {
        bined_groups.insert(grp);
      }
    }
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
  otpt << " Using ExactMGCancelator and " << n_samples << " points.\n";
  if (settings::use_virtual_collisions)
    otpt << "   Will use particles born after virtual collisions.\n";
  else
    otpt << "   Will NOT use particles born after virtual collisions.\n";
  Output::instance().write(otpt.str());

  return std::make_shared<ExactMGCancelator>(
      r_low, r_hi, shape, group_bins, settings::chi_matrix,
      settings::use_virtual_collisions, n_samples);
}
