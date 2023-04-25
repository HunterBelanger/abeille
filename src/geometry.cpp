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
#include <geometry/geometry.hpp>
#include <memory>
#include <utils/constants.hpp>
#include <utils/error.hpp>

namespace geometry {

//==========================================================================
// Global Variable Declrations
std::vector<std::shared_ptr<Surface>> surfaces;
std::vector<std::size_t> boundaries;
std::vector<std::shared_ptr<Cell>> cells;
std::shared_ptr<Universe> root_universe;
std::vector<std::shared_ptr<Universe>> universes;
std::vector<std::shared_ptr<Lattice>> lattices;

//==========================================================================
// Boundary struct constructor
Boundary::Boundary(double d, int index, BoundaryType bound)
    : distance{d}, surface_index(index), boundary_type{bound}, token(0) {}

//==========================================================================
// Function Definitions
Cell* get_cell(const Position& r, const Direction& u, int32_t on_surf) {
  // Ask root_universe for cell. If no cell is found, answer
  // will be a nullptr
  return root_universe->get_cell(r, u, on_surf);
}

Cell* get_cell(std::vector<GeoLilyPad>& stack, const Position& r,
               const Direction& u, int32_t on_surf) {
  // Ask root_universe for cell. If no cell is found, answer
  // will be a nullptr
  return root_universe->get_cell(stack, r, u, on_surf);
}

Boundary get_boundary(const Position& r, const Direction& u, int32_t on_surf) {
  // Define initial boundary which is first surface at INF
  double d_min = INF;
  BoundaryType btype = BoundaryType::Vacuum;
  int32_t closest_surface_indx = -1;

  // Go through all boundaries
  for (size_t i = 0; i < boundaries.size(); i++) {
    if (std::abs(on_surf) - 1 != static_cast<int>(boundaries[i])) {
      int32_t indx = static_cast<int32_t>(boundaries[i]);
      double dist =
          surfaces[static_cast<std::size_t>(indx)]->distance(r, u, false);
      if (dist < d_min) {
        closest_surface_indx = indx;
        d_min = dist;
      }
    }
  }

  if (closest_surface_indx >= 0) {
    return {
        d_min, closest_surface_indx,
        surfaces[static_cast<std::size_t>(closest_surface_indx)]->boundary()};
  }

  return {d_min, closest_surface_indx, btype};
}

int32_t id_to_token(int32_t id) {
  if (surface_id_to_indx.find(static_cast<uint32_t>(std::abs(id))) ==
      surface_id_to_indx.end())
    return 0;

  int32_t token = static_cast<int32_t>(
      surface_id_to_indx[static_cast<uint32_t>(std::abs(id))]);
  token += 1;

  return token;
}

void do_reflection(Particle& p, Boundary boundary) {
  // First, we should get the surface from the index
  if (boundary.surface_index < 0) {
    fatal_error("Surface index is less than zero in reflection.");
  }
  const std::shared_ptr<Surface>& surface =
      surfaces[static_cast<std::size_t>(boundary.surface_index)];

  // Get new Position object to temporarily contain the current position
  Position r_pre_refs = p.r();

  // In the event a reflection occurs imediately after another, we must
  // ensure that we use r to be the previous position, or the position before
  // the first reflection, as that is our previous true position which
  // resulted from a sampled collision
  if (p.is_reflected()) r_pre_refs = p.previous_r();
  Direction u = p.u();

  // Get position of particle on surface
  Position r_on_surf = p.r() + boundary.distance * u;

  Direction n = surface->norm(r_on_surf);  // Get norm of surface

  // Get new direction after reflection
  Vector new_dir = u - 2. * (u * n) * n;  // Calc new direction
  Direction u_new = Direction{new_dir.x(), new_dir.y(), new_dir.z()};

  // Calc new previous position before reflections
  double d = boundary.distance + (p.r() - r_pre_refs).norm();
  Position r_prev = r_on_surf - d * u_new;

  // Update particle MUST BE DONE IN THIS ORDER
  // as set_position will also change previous_r, so we must set that second
  p.set_position(r_on_surf);
  p.set_previous_r(r_prev);
  p.set_direction(u_new);
  p.set_reflected(true);
}

}  // namespace geometry
