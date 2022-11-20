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
#ifndef CELL_H
#define CELL_H

#include <limits>
#include <map>
#include <memory>
#include <vector>

#include <geometry/surfaces/surface.hpp>
#include <materials/material.hpp>

#include <yaml-cpp/yaml.h>

class Cell;

//============================================================================
// Operators for rpn cell definitions
enum OP : int32_t {
  L_PAR = std::numeric_limits<int32_t>::max(),
  R_PAR = std::numeric_limits<int32_t>::max() - 1,
  COMP = std::numeric_limits<int32_t>::max() - 2,
  INTR = std::numeric_limits<int32_t>::max() - 3,
  UNIN = std::numeric_limits<int32_t>::max() - 4
};

//============================================================================
// Externals from geometry.hpp
namespace geometry {
extern std::vector<std::shared_ptr<Surface>> surfaces;
extern std::vector<std::shared_ptr<Cell>> cells;
}  // namespace geometry

//============================================================================
// Eternals from parser.hpp
extern std::map<uint32_t, size_t> surface_id_to_indx;
extern std::map<uint32_t, size_t> cell_id_to_indx;

//============================================================================
// Cell Class
class Cell {
 public:
  Cell(std::vector<int32_t> i_rpn, std::shared_ptr<Material> material,
       uint32_t i_id, std::string i_name);
  ~Cell() = default;

  bool is_inside(const Position& r, const Direction& u, int32_t on_surf) const;

  std::pair<double, int32_t> distance_to_boundary(const Position& r,
                                                  const Direction& u,
                                                  int32_t on_surf) const;

  Material* material() { return material_raw_; }

  uint32_t id() const;

  std::string name() const;

  const std::vector<std::shared_ptr<Cell>> neighbors() const {
    return neighbors_;
  }

  void add_neighbor(std::shared_ptr<Cell> cell) { neighbors_.push_back(cell); }

 private:
  bool simple = true;
  std::vector<int32_t> rpn;  // Surface definition of cell
  uint32_t id_;
  std::string name_;

  bool is_inside_simple(const Position& r, const Direction& u,
                        int32_t on_surf) const;
  bool is_inside_complex(const Position& r, const Direction& u,
                         int32_t on_surf) const;

  std::shared_ptr<Material> material_;
  Material* material_raw_;

  std::vector<std::shared_ptr<Cell>> neighbors_{};

};  // Cell

//===========================================================================
// Non-Member Functions
std::vector<int32_t> infix_to_rpn(const std::vector<int32_t>& infix);

void make_cell(YAML::Node cell_node);

#endif
