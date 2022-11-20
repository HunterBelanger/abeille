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
#ifndef EXACT_MG_CANCELATOR_H
#define EXACT_MG_CANCELATOR_H

#include <array>
#include <cstdint>
#include <materials/material.hpp>
#include <materials/mg_nuclide.hpp>
#include <optional>
#include <simulation/cancelator.hpp>
#include <unordered_map>
#include <vector>

class ExactMGCancelator : public Cancelator {
 public:
  ExactMGCancelator(const Position& r_low, const Position& r_hi,
                    const std::array<std::size_t, 4>& shape,
                    const std::vector<std::vector<std::size_t>>& group_bins,
                    bool chi_matrix, bool use_virtual_collisions,
                    uint32_t n_samples);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation(pcg32&) override final;
  std::vector<BankedParticle> get_new_particles(pcg32& rng) override final;
  void clear() override final;

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

  // Key which represents a unique cancellation bin for a
  // given position and energy group.
  struct Key {
    Key(std::size_t i, std::size_t j, std::size_t k, std::size_t e)
        : i(i), j(j), k(k), e(e) {}

    std::size_t i, j, k, e;

    std::size_t hash_key() const {
      return e + shape[3] * (k + shape[2] * (j + shape[1] * i));
    }

    bool operator==(const Key& other) const {
      return ((i == other.i) && (j == other.j) && (k == other.k) &&
              (e == other.e));
    }

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

  struct KeyHash {
    std::size_t operator()(const Key& key) const { return key.hash_key(); }
  };

  //==========================================================================
  // Data Members

  // All cancellation bins, organized first by Key has, then by material.
  std::unordered_map<Key, std::unordered_map<Material*, CancelBin>, KeyHash>
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
                                          pcg32& rng) const;

  // Computes <f> and <1/f> for all particles in all bins.
  void compute_averages(const Key& key, Material* mat, MGNuclide* nuclide,
                        CancelBin& bin);

  void cancel_bin(CancelBin& bin, MGNuclide* nuclide, bool first_wgt);
};

std::shared_ptr<ExactMGCancelator> make_exact_mg_cancelator(
    const YAML::Node& node);

#endif
