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
#ifndef MG_NUCLIDE_H
#define MG_NUCLIDE_H

#include <yaml-cpp/yaml.h>

#include <cstdint>
#include <materials/mg_angle_distribution.hpp>
#include <materials/nuclide.hpp>
#include <vector>

class MGNuclide : public Nuclide {
 public:
  MGNuclide(const std::vector<double>& speeds, const std::vector<double>& Et,
            const std::vector<double>& Ea, const std::vector<double>& Ef,
            const std::vector<double>& nu_prmpt,
            const std::vector<double>& nu_dlyd,
            const std::vector<std::vector<double>>& chi,
            const std::vector<std::vector<double>>& Es,
            const std::vector<std::vector<double>>& yield,
            const std::vector<std::vector<MGAngleDistribution>>& angle,
            const std::vector<double>& P_dlyd_grp,
            const std::vector<double>& decay_cnsts);

  bool fissile() const override final;
  bool has_urr() const override final;
  double total_xs(double E_in, std::size_t i) const override final;
  double disappearance_xs(double E_in, std::size_t i) const override final;
  double fission_xs(double E_in, std::size_t i) const override final;
  double nu_total(double E_in, std::size_t i) const override final;
  double nu_prompt(double E_in, std::size_t i) const override final;
  double nu_delayed(double E_in, std::size_t i) const override final;
  double reaction_xs(uint32_t mt, double E_in, size_t i) const override final;
  double elastic_xs(double E_in, std::size_t i) const override final;
  std::size_t energy_grid_index(double E) const override final;
  MicroXSs get_micro_xs(double E, std::optional<double> urr_rand = std::nullopt) const override final;

  std::size_t num_delayed_groups() const override final;
  double delayed_group_constant(std::size_t g) const override final;
  double delayed_group_probability(std::size_t g,
                                   double E) const override final;

  ScatterInfo sample_scatter(double Ein, const Direction& u, const MicroXSs& micro_xs,
                             pcg32& rng) const override final;
  ScatterInfo sample_scatter_mt(uint32_t mt, double Ein, const Direction& u,
                                std::size_t i, pcg32& rng) const override final;
  FissionInfo sample_fission(double Ein, const Direction& u, std::size_t i,
                             double Pdelayed, pcg32& rng) const override final;
  FissionInfo sample_prompt_fission(double Ein, const Direction& u,
                                    std::size_t i,
                                    pcg32& rng) const override final;
  FissionInfo sample_delayed_fission(double Ein, const Direction& u,
                                     std::size_t g,
                                     pcg32& rng) const override final;

  double max_energy() const override final;
  double min_energy() const override final;
  double speed(double E, std::size_t i) const override final;

  uint32_t zaid() const override final;

  // Methods unique to MGNuclide, for exact MG cancellation.
  const std::vector<double>& group_speeds() const { return group_speeds_; }
  const std::vector<double>& Et() const { return Et_; }
  const std::vector<double>& Ea() const { return Ea_; }
  const std::vector<double>& Ef() const { return Ef_; }
  const std::vector<double>& Es() const { return Es_; }
  const std::vector<double>& nu_prmpt() const { return nu_prmpt_; }
  const std::vector<double>& nu_dlyd() const { return nu_delyd_; }
  const std::vector<std::vector<double>>& chi() const { return chi_; }
  const std::vector<std::vector<double>>& Ps() const { return Ps_; }
  const std::vector<std::vector<double>>& mult() const { return mult_; }
  const std::vector<std::vector<MGAngleDistribution>>& angles() const {
    return angle_dists_;
  }
  const std::vector<double>& P_dlyd_grp() const { return P_delayed_group; }
  const std::vector<double>& dlyd_grp_decay_const() const {
    return delayed_group_decay_constants;
  }

 private:
  std::vector<double> group_speeds_;
  std::vector<double> Et_;
  std::vector<double> Ea_;
  std::vector<double> Ef_;
  std::vector<double> Es_;
  std::vector<double> nu_prmpt_;
  std::vector<double> nu_delyd_;
  std::vector<std::vector<double>> chi_;
  std::vector<std::vector<double>> Ps_;
  std::vector<std::vector<double>> mult_;
  std::vector<std::vector<MGAngleDistribution>> angle_dists_;
  std::vector<double> P_delayed_group;
  std::vector<double> delayed_group_decay_constants;
  bool fissile_ = false;

  void verify() const;
  void make_scatter_xs();
  void normalize_chi();
  void check_sizes() const;
  void check_xs() const;
  void check_fission_data() const;
  void check_dealyed_data() const;
  void check_fissile();
};

std::shared_ptr<MGNuclide> make_mg_nuclide(const YAML::Node& mat, uint32_t id);

#endif
