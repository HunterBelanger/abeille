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
#ifndef SLICE_PLOT_H
#define SLICE_PLOT_H

#include <geometry/cell.hpp>
#include <map>
#include <mutex>
#include <plotting/pixel.hpp>
#include <string>
#include <utils/direction.hpp>
#include <utils/position.hpp>
#include <utils/rng.hpp>
#include <vector>

namespace plotter {
// Maps for colors from plotter.hpp
extern std::map<uint32_t, Pixel> cell_id_to_color;
extern std::map<uint32_t, Pixel> material_id_to_color;

class SlicePlot {
 public:
  enum Basis { XY, XZ, YZ };
  enum ColorBy { Cell, Material };

  SlicePlot(std::string fname, uint64_t pwidth, uint64_t pheight, double width,
            double height, Position origin, Basis basis, ColorBy colorby,
            Pixel background = Pixel(255, 255, 255));
  ~SlicePlot() = default;

  void generate_plot();
  void write() const;

 private:
  std::string file_name_;
  uint64_t plot_width_, plot_height_;  // Both are in number of pixels
  Basis basis_;
  ColorBy colorby_;
  Position origin_;
  Pixel background_;
  std::vector<Pixel> image_matrix;  // Row-Major order
  double width_, height_;
  pcg32 rng;
  std::mutex create_color_mutex;

  Position get_pixel_position(uint64_t i, uint64_t j) const;
  Position get_start_position(uint64_t i) const;
  Direction get_tracking_direction() const;
  Pixel get_color(::Cell* cell);
  Pixel get_random_color();

};  // SlicePlot

}  // namespace plotter

#endif
