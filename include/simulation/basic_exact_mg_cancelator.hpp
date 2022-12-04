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
#ifndef BASIC_EXACT_MG_CANCELATOR_H
#define BASIC_EXACT_MG_CANCELATOR_H

#include <functional>
#include <materials/material_helper.hpp>
#include <memory>
#include <optional>
#include <simulation/cancelator.hpp>
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
