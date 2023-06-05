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
#ifndef ENTROPY_H
#define ENTROPY_H

#include <array>
#include <cmath>
#include <cstdint>
#include <utils/position.hpp>
#include <vector>

class Entropy {
 public:
  enum Sign { Positive, Negative, Total };

 public:
  Entropy(Position low, Position up, std::array<uint32_t, 3> shp, Sign sgn)
      : lower_corner{low},
        upper_corner{up},
        shape{shp},
        nbins{0},
        dx{},
        dy{},
        dz{},
        bins{},
        total_weight{},
        sign{sgn} {
    nbins = shape[0] * shape[1] * shape[2];
    bins.resize(nbins, 0.);
    total_weight = 0.;
    dx = (upper_corner.x() - lower_corner.x()) / static_cast<double>(shape[0]);
    dy = (upper_corner.y() - lower_corner.y()) / static_cast<double>(shape[1]);
    dz = (upper_corner.z() - lower_corner.z()) / static_cast<double>(shape[2]);
  }
  ~Entropy() = default;

  void add_point(const Position& r, const double& w);

  void synchronize_entropy_across_nodes();

  double calculate_entropy() const;

  double calculate_empty_fraction() const;

  void zero();

  double total_wgt() const { return total_weight; }

 private:
  Position lower_corner;
  Position upper_corner;
  std::array<std::uint32_t, 3> shape;
  std::uint32_t nbins;
  double dx, dy, dz;
  std::vector<double> bins;
  double total_weight;
  Sign sign;

};  // Etnropy

#endif  // MG_ENTROPY_H
