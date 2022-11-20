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
#ifndef HEX_LATTICE_H
#define HEX_LATTICE_H

#include <array>
#include <geometry/lattice.hpp>
#include <utils/constants.hpp>

class HexLattice : public Lattice {
 public:
  enum Top { Pointy, Flat };

  HexLattice(uint32_t nrings, uint32_t nz, double p, double pz, double x,
             double y, double z, Top t, uint32_t i_id, std::string i_name);
  ~HexLattice() = default;

  bool is_inside(Position r, Direction u) const override;

  std::array<int32_t, 3> get_tile(Position p, Direction u) const override;

  Cell* get_cell(Position r, Direction u, int32_t on_surf) const override;

  Cell* get_cell(std::vector<GeoLilyPad>& stack, Position r, Direction u,
                 int32_t on_surf) const override;

  double distance_to_tile_boundary(Position r_local, Direction u,
                                   std::array<int32_t, 3> tile) const override;

  void set_elements(std::vector<int32_t> univs) override;

 private:
  double pitch_, pitch_z_;
  double X_o, Y_o, Z_o;  // Origin of center hex
  const double cos_pi_6 = std::cos(PI / 6.0);
  const double sin_pi_6 = std::sin(PI / 6.0);
  const double cos_pi_3 = std::cos(PI / 3.0);
  const double sin_pi_3 = std::sin(PI / 3.0);
  uint32_t Nrings, Nz, Nhex;
  uint32_t width, mid_qr;
  Top top_;

  std::array<int32_t, 2> get_nearest_hex(Position p) const;

  Position get_hex_center(std::array<int32_t, 2> qr) const;

  Position tile_center(int q, int r, int nz) const;

  uint32_t get_ring(std::array<int32_t, 2> qr) const;

  double distance_to_tile_boundary_pointy(Position r_tile, Direction u) const;

  double distance_to_tile_boundary_flat(Position r_tile, Direction u) const;

  double distance_to_line(Position r, Direction u, double x1, double y1,
                          double x2, double y2) const;

  size_t linear_index(std::array<int32_t, 2> qr, uint32_t z) const;

};  // HexLattice

void make_hex_lattice(YAML::Node latt_node, YAML::Node input);

#endif
