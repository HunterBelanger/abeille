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
#ifndef CELL_SEARCH_MESH_H
#define CELL_SEARCH_MESH_H

#include <yaml-cpp/yaml.h>

#include <array>
#include <geometry/cell.hpp>

class CellSearchMesh {
 public:
  CellSearchMesh(
      Position low, Position hi, const std::array<uint32_t, 3>& shape,
      const std::vector<std::vector<std::shared_ptr<Cell>>>& elements);

  bool is_inside(const Position& r, const Direction& u) const;
  std::size_t index(const Position& r, const Direction& u) const;
  const std::vector<std::shared_ptr<Cell>>& index_cells(std::size_t i) const;
  Cell* find_cell(const Position& r, const Direction& u, int32_t on_surf) const;
  std::size_t size() const { return elements_.size(); }

 private:
  std::array<uint32_t, 3> shape_;
  std::array<double, 3> pitch_;
  Position r_low, r_hi;
  std::vector<std::vector<std::shared_ptr<Cell>>> elements_;

  std::array<int32_t, 3> get_tile(const Position& r, const Direction& u) const;
};

std::shared_ptr<CellSearchMesh> make_cell_search_mesh(const YAML::Node& node);

#endif
