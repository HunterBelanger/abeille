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
#ifndef NUCLIDE_H
#define NUCLIDE_H

#include <cstdint>
#include <simulation/noise_source.hpp>
#include <simulation/particle.hpp>

#include <map>
#include <optional>
#include <unordered_set>

class Nuclide;
extern std::map<uint32_t, std::shared_ptr<Nuclide>> nuclides;
extern std::unordered_set<uint32_t> zaids_with_urr;

struct MicroXSs {
  double total = 0.;
  double fission = 0.;
  double absorption = 0.;
  double elastic = 0.;
  double inelastic = 0.; // Not actually inelastic, just all scattering that isn't elastic
  double nu_total = 0.;
  double nu_delayed = 0.;
  double energy = 0.;
  double concentration = 0.;
  double noise_copy = 0.;
  std::size_t energy_index = 0;
  bool urr = false;
};

struct ScatterInfo {
  double yield = 1.;
  double energy = 0.;
  Direction direction = Direction(1., 0., 0.);
  uint32_t mt = 0;
  bool thermal = false;
};

struct FissionInfo {
  double energy = 0.;
  Direction direction = Direction(1., 0., 0.);
  bool delayed = false;
  double precursor_decay_constant = 0.;
  uint32_t delayed_family = 0;
};

// This is a general nuclide interface, to permit use with
// the same transport and simulation classes, regardless of
// the energy mode being used.
class Nuclide {
 public:
  Nuclide() : id_(id_counter++) {}
  virtual ~Nuclide() = default;

  virtual bool fissile() const = 0;
  virtual bool has_urr() const = 0;
  virtual double total_xs(double E_in, std::size_t i) const = 0;
  virtual double disappearance_xs(double E_in, std::size_t i) const = 0;
  virtual double fission_xs(double E_in, std::size_t i) const = 0;
  virtual double nu_total(double E_in, std::size_t i) const = 0;
  virtual double nu_prompt(double E_in, std::size_t i) const = 0;
  virtual double nu_delayed(double E_in, std::size_t i) const = 0;
  virtual double reaction_xs(uint32_t mt, double E_in, size_t i) const = 0;
  virtual double elastic_xs(double E_in, std::size_t i) const = 0;
  virtual std::size_t energy_grid_index(double E) const = 0;
  virtual MicroXSs get_micro_xs(double E, std::optional<double> urr_rand = std::nullopt) const = 0;

  virtual std::size_t num_delayed_groups() const = 0;
  virtual double delayed_group_constant(std::size_t g) const = 0;
  virtual double delayed_group_probability(std::size_t g, double E) const = 0;

  virtual ScatterInfo sample_scatter(double Ein, const Direction& u,
                                     const MicroXSs& micro_xs, pcg32& rng) const = 0;
  virtual ScatterInfo sample_scatter_mt(uint32_t mt, double Ein,
                                        const Direction& u, std::size_t i,
                                        pcg32& rng) const = 0;
  virtual FissionInfo sample_fission(double Ein, const Direction& u,
                                     std::size_t i, double Pdelayed,
                                     pcg32& rng) const = 0;
  virtual FissionInfo sample_prompt_fission(double Ein, const Direction& u,
                                            std::size_t i,
                                            pcg32& rng) const = 0;
  // Samples an energy and direction from the delayed spectrum of delayed family
  // g.
  virtual FissionInfo sample_delayed_fission(double Ein, const Direction& u,
                                             std::size_t g,
                                             pcg32& rng) const = 0;

  virtual double max_energy() const = 0;
  virtual double min_energy() const = 0;
  virtual double speed(double E, std::size_t i) const = 0;

  virtual uint32_t zaid() const = 0;

  uint32_t id() const { return id_; }

 private:
  static uint32_t id_counter;
  uint32_t id_;
};

#endif
