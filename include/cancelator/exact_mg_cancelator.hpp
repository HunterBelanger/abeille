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
#ifndef EXACT_MG_CANCELATOR_H
#define EXACT_MG_CANCELATOR_H

#include <cancelator/cancelator.hpp>
#include <materials/material.hpp>
#include <materials/mg_nuclide.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

class ExactMGCancelator : public Cancelator {
 public:
  ExactMGCancelator(const Position& r_low, const Position& r_hi,
                    const std::array<uint64_t, 4>& shape,
                    const std::vector<std::vector<std::size_t>>& group_bins,
                    bool chi_matrix, bool use_virtual_collisions,
                    uint32_t n_samples);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation() override final;
  std::vector<BankedParticle> get_new_particles(RNG& rng) override final;
  void clear() override final;
  void check_particle_mover_compatibility(
      const std::shared_ptr<IParticleMover>& pmover) const override final;

  void write_output_info(H5::Group& grp) const override final;

  // Key which represents a unique cancellation bin for a
  // given position and energy group.
  // Key is public for MPI use
  struct Key {
    Key(uint64_t i, uint64_t j, uint64_t k, uint64_t e)
        : i(i), j(j), k(k), e(e) {}

    Key() : i(0), j(0), k(0), e(0) {}
    uint64_t i, j, k, e;

    std::size_t hash_key() const {
      return e + shape[3] * (k + shape[2] * (j + shape[1] * i));
    }

    bool operator==(const Key& other) const {
      return ((i == other.i) && (j == other.j) && (k == other.k) &&
              (e == other.e));
    }

    auto operator<=>(const Key&) const = default;

    // Contains the shape of the cancellation region mesh.
    // shape[0] Number of regions in x
    // shape[1] Number of regions in y
    // shape[2] Number of regions in z
    // shape[3] Number of regions in energy
    static std::array<std::size_t, 4> shape;

    // The width of each region in x, y, and z
    static std::array<double, 3> pitch;

    // Contains the groupings for the energy bins.
    // goup_bin.size() == shape[3] should ALWAYS be true !
    // group_bin[i] is then a vector containing the energy
    // group indices which belong to the ith energy bin.
    static std::vector<std::vector<std::size_t>> group_bins;

    static Position r_low, r_hi;
  };

 private:
  //==========================================================================
  // Objects
  struct CancelBin {
    struct Averages {
      double f = 0.;
      double f_inv = 0.;
    };

    double uniform_wgt = 0.;
    double uniform_wgt2 = 0.;
    double W = 0.;
    double W2 = 0.;
    double sum_c = 0.;
    double sum_c_wgt = 0.;
    double sum_c_wgt2 = 0.;
    bool can_cancel = true;
    std::vector<BankedParticle*> particles;
    std::vector<Averages> averages;
  };

  struct KeyHash {
    std::size_t operator()(const Key& key) const { return key.hash_key(); }
  };

  //==========================================================================
  // Data Members

  // All cancellation bins, organized first by Key has, then by material.
  std::unordered_map<Key, std::unordered_map<uint32_t, CancelBin>, KeyHash>
      bins;

  // False if we only have MG materials with a chi vector.
  const bool CHI_MATRIX;
  // Number of samples to use when computing <f> and <1/f>.
  const uint32_t N_SAMPLES;
  // Max number of attempts to sample a position within the region / material
  // before we give up and say we can't cancel that particle.
  const uint32_t N_MAX_POS{100};

  // If false, we don't cancel particles who's parents had a virtual collision
  // just before, as we cannot guarentee the exactness of the method.
  const bool USE_VIRTUAL_COLLISIONS;

  //==========================================================================
  // Private Methods

  std::optional<Key> get_key(const Position& r, std::size_t g);

  std::vector<std::pair<ExactMGCancelator::Key, uint32_t>> sync_keys();

  // Get's a pointer to the material at r
  Material* get_material(const Position& r) const;

  // Calculates the transition kernel from P1 to P4, where we assume that
  // the fission angular distribution at P4 is perfectly isotropic.
  double get_f(const Position& r1, const Direction& u1, std::size_t g1,
               std::size_t g3, const Position& r4, std::size_t g4, double Esmp,
               MGNuclide* nuclide) const;

  double get_beta(const CancelBin& bin, std::size_t i, bool wgt_1) const;

  std::optional<std::pair<Position, std::size_t>> sample_point(
      const Key& key, Material* mat, unsigned long long& i) const;

  std::optional<Position> sample_position(const Key& key, Material* mat,
                                          RNG& rng) const;

  // Computes <f> and <1/f> for all particles in all bins.
  void compute_averages(const Key& key, Material* mat, MGNuclide* nuclide,
                        CancelBin& bin);

  void cancel_bin(CancelBin& bin, MGNuclide* nuclide);
};

std::shared_ptr<ExactMGCancelator> make_exact_mg_cancelator(
    const YAML::Node& node);

#endif
