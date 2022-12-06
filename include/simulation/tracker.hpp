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
#ifndef TRACKER_H
#define TRACKER_H

#include <geometry/cell.hpp>
#include <geometry/geo_lily_pad.hpp>
#include <geometry/geometry.hpp>
#include <materials/material.hpp>
#include <simulation/particle.hpp>
#include <sstream>
#include <stdexcept>
#include <utils/error.hpp>
#include <utils/parser.hpp>

class Tracker {
 public:
  Tracker(Position i_r, Direction i_u, int32_t token = 0)
      : r_(i_r), u_(i_u), tree(), surface_token_(token) {
    // Allocate enough space for at least 10 levels
    tree.reserve(10);

    current_cell = geometry::get_cell(tree, r_, u_, surface_token_);
    if (current_cell) current_mat = current_cell->material();
  };

  ~Tracker() = default;

  Position r() const { return r_; }
  Direction u() const { return u_; }

  void set_r(Position r) {
    r_ = r;
    surface_token_ = 0;
  }
  void set_u(Direction u) { u_ = u; }

  void restart_get_current() {
    tree.clear();
    current_cell = geometry::get_cell(tree, r_, u_, surface_token_);
    if (current_cell)
      current_mat = current_cell->material();
    else
      current_mat = nullptr;
  }

  void move(double d) {
    r_ = r_ + d * u_;

    // Go update positions inside tree
    for (auto& leaf : tree) {
      leaf.r_local = leaf.r_local + d * u_;
    }

    surface_token_ = 0;
  }

  int32_t surface_token() const { return surface_token_; }
  void set_surface_token(int32_t st) { surface_token_ = st; }

  Material* material() { return current_mat; }
  Cell* cell() { return current_cell; }

  geometry::Boundary boundary() const {
    return geometry::get_boundary(r_, u_, surface_token_);
  }

  // geometry::Boundary __attribute__ ((noinline)) get_nearest_boundary() const
  // {
  geometry::Boundary get_nearest_boundary() const {
    auto bound = this->boundary();

    if (!this->is_lost()) {
      double dist = bound.distance;
      BoundaryType btype = bound.boundary_type;
      int surface_index = bound.surface_index;
      int32_t token = bound.token;

      // Go up the entire tree
      for (const auto& pad : tree) {
        if (pad.type == GeoLilyPad::PadType::Lattice) {
          auto lat_id = lattice_id_to_indx[pad.id];
          double d = geometry::lattices[lat_id]->distance_to_tile_boundary(
              pad.r_local, u_, pad.tile);
          if (d < dist && std::abs(d - dist) > 100 * SURFACE_COINCIDENT) {
            dist = d;
            btype = BoundaryType::Normal;
            surface_index = -1;
            token = 0;
          }
        } else if (pad.type == GeoLilyPad::PadType::Cell) {
          auto cell_id = cell_id_to_indx[pad.id];
          auto d_t = geometry::cells[cell_id]->distance_to_boundary(
              pad.r_local, u_, surface_token_);
          if (d_t.first < dist &&
              std::abs(d_t.first - dist) > 100 * SURFACE_COINCIDENT) {
            dist = d_t.first;
            token = std::abs(d_t.second);

            if (token) {
              surface_index = token - 1;
            } else {
              surface_index = -1;
            }

            if (surface_index >= 0)
              btype =
                  geometry::surfaces[static_cast<std::size_t>(surface_index)]
                      ->boundary();
            else
              btype = BoundaryType::Normal;

            if (surface_index >= 0 &&
                geometry::surfaces[static_cast<std::size_t>(surface_index)]
                        ->sign(pad.r_local, u_) < 0)
              token *= -1;
          }
        }
      }

      geometry::Boundary ret_bound(dist, surface_index, btype);
      ret_bound.token = token;

      // Distance to surface, and token, with sign indicating the positions
      // current orientation to the surface.
      return ret_bound;
    } else {
      return geometry::get_boundary(r_, u_, surface_token_);
    }
  }

  void cross_surface(geometry::Boundary d_t) {
    this->move(d_t.distance);
    // Use negative as we have changed orientation !
    surface_token_ = -d_t.token;
  }

  bool is_lost() const { return !current_cell; }

  /*void get_current() {
    current_cell = geometry::get_cell(r_, u_, surface_token_);
    current_mat = current_cell->material();
  }*/

  void get_current() {
    // Iterator pointing to the highest leaf in the tree
    // which is not true (either cell or lattice)
    auto first_bad = tree.end();

    if (!check_tree()) {
      // Start from scratch
      restart_get_current();

      if (!check_tree()) {
        std::stringstream mssg;
        mssg << " BAD POSITIONS!\n r_ = " << r_
             << "\n rl = " << tree.front().r_local << "\n";
        fatal_error(mssg.str(), __FILE__, __LINE__);
      }
    }

    // Go back through tree, and see where we are no-longer inside
    for (auto it = tree.begin(); it != tree.end(); it++) {
      if (it->type == GeoLilyPad::PadType::Cell) {
        auto cell_indx = cell_id_to_indx[it->id];
        const auto& cell = geometry::cells[cell_indx];
        if (!cell->is_inside(it->r_local, u_, surface_token_)) {
          first_bad = it;
          break;
        }
      } else if (it->type == GeoLilyPad::PadType::Lattice) {
        auto lat_indx = lattice_id_to_indx[it->id];
        const auto& lat = geometry::lattices[lat_indx];
        auto tile = lat->get_tile(it->r_local, u_);
        // Check if tile has changed
        if (it->tile[0] != tile[0] || it->tile[1] != tile[1] ||
            it->tile[2] != tile[2]) {
          first_bad = it;
          break;
        }
      }
    }

    // Only need to research if last_bad != tree.rend(), otherwise the
    // cell and material have not changed.
    if (first_bad != tree.end()) {
      // Get rid of bad tree elements. Get index of last bad, which is
      // the size of the number of good elements.
      auto size =
          static_cast<std::size_t>(std::distance(tree.begin(), first_bad));
      tree.resize(size);

      // Now start at the last element, and get the new position
      auto uni_indx = universe_id_to_indx[tree.back().id];
      const auto& uni = geometry::universes[uni_indx];
      Position r_local = tree.back().r_local;
      tree.pop_back();
      current_cell = uni->get_cell(tree, r_local, u_, surface_token_);
      if (!current_cell) restart_get_current();
      if (current_cell) current_mat = current_cell->material();
    }
  }

  bool check_tree() const {
    return r_.x() == tree.front().r_local.x() &&
           r_.y() == tree.front().r_local.y() &&
           r_.z() == tree.front().r_local.z();
  }

  void do_reflection(Particle& p, geometry::Boundary boundary) {
    // Get the surface first
    if (boundary.surface_index < 0) {
      fatal_error("Bad surface index in Tracker::do_reflection", __FILE__,
                  __LINE__);
    }
    const std::shared_ptr<Surface>& surface =
        geometry::surfaces[static_cast<std::size_t>(boundary.surface_index)];

    int32_t token = geometry::id_to_token(static_cast<int32_t>(surface->id()));
    if (surface->sign(p.r(), p.u()) < 0) token *= -1;
    this->set_surface_token(token);

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

    this->set_r(p.r());
    this->set_u(p.u());
    this->restart_get_current();
  }

 private:
  Position r_;
  Direction u_;
  std::vector<GeoLilyPad> tree;
  Material* current_mat = nullptr;
  Cell* current_cell = nullptr;
  // Token for current surface which particle is on. The index for the
  // surface in the surface vector is std::abs(token) - 1. Token is NOT
  // the surface id !!
  int32_t surface_token_ = 0;

};  // Tracker

#endif  // MG_TRACKER_H
