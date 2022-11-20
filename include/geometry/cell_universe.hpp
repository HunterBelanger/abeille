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
#ifndef CELL_UNIVERSE_H
#define CELL_UNIVERSE_H

#include <yaml-cpp/yaml.h>

#include <geometry/universe.hpp>
#include <map>

//===========================================================================
// Externals from geometry
namespace geometry {
extern std::vector<std::shared_ptr<Cell>> cells;
extern std::vector<std::shared_ptr<Universe>> universes;
}  // namespace geometry

//===========================================================================
// Externals from parser
extern std::map<uint32_t, size_t> cell_id_to_indx;
extern std::map<uint32_t, size_t> universe_id_to_indx;

//===========================================================================
// CellUniverse class
class CellUniverse : public Universe {
 public:
  CellUniverse(std::vector<uint32_t> i_ind, uint32_t i_id, std::string i_name);
  ~CellUniverse() = default;

  Cell* get_cell(Position r, Direction u, int32_t on_surf) const override;

  Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                 int32_t on_surf) const override;

 private:
  std::vector<uint32_t> cell_indicies;

};  // CellUniverse

//===========================================================================
// Non-Member functions
void make_cell_universe(YAML::Node uni_node);

#endif
