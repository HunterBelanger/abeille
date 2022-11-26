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
#ifndef PIXEL_H
#define PIXEL_H

#include <cstdint>
#include <iostream>

namespace plotter {

//======================================================================
// Pixel class to contain info for a pixel or color definition
class Pixel {
 public:
  Pixel(uint8_t i_r = 0xff, uint8_t i_g = 0xff, uint8_t i_b = 0xff): r(i_r), g(i_g), b(i_b) {}
  ~Pixel() = default;

  uint8_t R() const { return r; }
  uint8_t G() const { return g; }
  uint8_t B() const { return b;}

  void set_R(uint8_t i_r) { r = i_r; }
  void set_G(uint8_t i_g) { g = i_g; }
  void set_B(uint8_t i_b) { b = i_b; }

 private:
  uint8_t r, g, b;
};  // Pixel

//========================================================================
// Non-Member Function Declarations
inline std::ostream& operator<<(std::ostream& out, const Pixel& p) {
  out << p.R() << p.G() << p.B();
  return out;
}

}  // namespace plotter

#endif
