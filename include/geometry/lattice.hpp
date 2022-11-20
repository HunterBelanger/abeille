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
#ifndef LATTICE_H
#define LATTICE_H

#include <yaml-cpp/yaml.h>

#include <geometry/universe.hpp>
#include <map>

class Lattice;

//===========================================================================
// Externals from geometry
namespace geometry {
extern std::vector<std::shared_ptr<Lattice>> lattices;
extern std::vector<std::shared_ptr<Universe>> universes;
}  // namespace geometry

//===========================================================================
// Externals from parser
extern std::map<uint32_t, size_t> universe_id_to_indx;
extern std::map<uint32_t, size_t> lattice_id_to_indx;
extern void find_universe(YAML::Node input, uint32_t id);

class Lattice {
 public:
  Lattice(uint32_t i_id, std::string i_name);
  virtual ~Lattice() = default;

  // Returns ture if position is inside a true lattice element
  // (not nullptr element), and false if it is outside or in a
  // dummy lattice element
  virtual bool is_inside(Position r, Direction u) const = 0;

  virtual std::array<int32_t, 3> get_tile(Position r, Direction u) const = 0;

  // Finds lattice element containing given position, transforms
  // coordinated to that elements frame, then asks that universe
  // for the cell of the local coordiante given.
  virtual Cell* get_cell(Position r, Direction u, int32_t on_surf) const = 0;

  virtual Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r,
                         Direction u, int32_t on_surf) const = 0;

  // Given the position in the frame of the lattice (NOT THE FRAME OF THE
  // TILE!), the distance to the edge of the provided tile is returned.
  virtual double distance_to_tile_boundary(
      Position r_local, Direction u, std::array<int32_t, 3> tile) const = 0;

  virtual void set_elements(std::vector<int32_t> univs) = 0;

  void set_outisde_universe(int32_t univ);

  uint32_t id() const;

  std::string name() const;

 protected:
  // Any lattice element that is negative gets directed to the
  // outer_universe. If outer_universe_index is also negative,
  // the particle is considered to be lost.
  std::vector<int32_t> lattice_universes;
  int32_t outer_universe_index;
  uint32_t id_;
  std::string name_;

};  // Lattice

#endif
