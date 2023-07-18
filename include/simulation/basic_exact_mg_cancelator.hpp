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
#ifndef BASIC_EXACT_MG_CANCELATOR_H
#define BASIC_EXACT_MG_CANCELATOR_H

#include <materials/material_helper.hpp>
#include <simulation/cancelator.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <unordered_map>

class BasicExactMGCancelator : public Cancelator {
 public:
  enum class BetaMode { Zero, Minimum, OptAverageF, OptAverageGain };

  BasicExactMGCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                         uint32_t Nz, BetaMode beta, bool sobol, uint32_t nsmp);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation(pcg32& rng) override final;
  std::vector<BankedParticle> get_new_particles(pcg32& rng) override final;
  void clear() override final;

 private:
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
    uint64_t rng_seed_advance = 0;
    bool can_cancel = true;
    std::vector<BankedParticle*> particles;
    std::vector<Averages> averages;
  };

  struct Key {
    int i, j, k;

    bool operator==(const Key& other) const {
      return ((i == other.i) && (j == other.j) && (k == other.k));
    }
  };

  class KeyHash {
   public:
    static std::array<uint32_t, 3> shape;
    std::size_t operator()(const Key& key) const {
      int int_key =
          key.k + static_cast<int>(this->shape[2]) *
                      (key.j + static_cast<int>(this->shape[1]) * key.i);
      return std::hash<int>()(int_key);
    }
  };

  const Position r_low, r_hi;
  KeyHash hash_fn;
  const double dx, dy, dz;
  const BetaMode beta_mode;
  const bool use_sobol;
  std::unordered_map<Key, std::unordered_map<Material*, CancelBin>, KeyHash>
      bins;
  const uint32_t N_SAMPLES;
  const uint32_t N_MAX_POS = 100;  // Max number of position samples

  //==========================================================================
  // Private Helper Methods

  void get_averages(const Key& key, Material* mat, CancelBin& bin, pcg32& rng);

  void get_averages_sobol(const Key& key, Material* mat, CancelBin& bin);

  Material* get_material(const Position& r) const;

  std::optional<Position> sample_position(const Key& key, Material* mat,
                                          pcg32& rng) const;

  std::optional<Position> sample_position_sobol(const Key& key, Material* mat,
                                                unsigned long long& i) const;

  double get_f(const Position& r, const Position& r_parent, double Esmp) const;

  double get_min_f(const Key& key, const Position& r_parent, double Esmp) const;

  double get_beta(const Key& key, const CancelBin& bin, std::size_t i,
                  const Position& r_parent, double Esmp, double wgt,
                  bool first_wgt) const;

  void cancel_bin(const Key& key, Material* mat, CancelBin& bin,
                  bool first_wgt);
};

std::shared_ptr<BasicExactMGCancelator> make_basic_exact_mg_cancelator(
    const YAML::Node& node);

#endif
