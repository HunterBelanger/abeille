/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef TRACKER_H
#define TRACKER_H

#include <geometry/boundary.hpp>
#include <geometry/cell.hpp>
#include <geometry/geo_lily_pad.hpp>
#include <geometry/geometry.hpp>
#include <materials/material.hpp>
#include <simulation/particle.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/parser.hpp>

#include <sstream>
#include <stdexcept>

class Tracker {
 public:
  Tracker(Position i_r, Direction i_u, int32_t token = 0)
      : r_(i_r), u_(i_u), tree(), surface_token_(token) {
    // Allocate enough space for at least 10 levels
    tree.reserve(10);

    current_cell = geometry::get_cell(tree, r_, u_, surface_token_);
    if (current_cell) current_mat = current_cell.cell->material();
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
    if (current_cell) {
      if (current_cell.cell->fill() == Cell::Fill::Universe) {
        fatal_error("Did not find a cell with a material.");
      }
      current_mat = current_cell.cell->material();
    } else {
      current_mat = nullptr;
    }
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
  Cell* cell() { return current_cell.cell; }
  uint32_t cell_instance() const { return current_cell.instance; }

  Boundary get_boundary_condition() const {
    if (!this->is_lost()) {
      double dist = INF;
      BoundaryType btype = BoundaryType::Vacuum;
      int surface_index = -1;
      int32_t token = 0;

      // Go up the entire tree
      for (const auto& pad : tree) {
        if (pad.type == GeoLilyPad::PadType::Cell) {
          auto cell_id = cell_id_to_indx[pad.id];
          Cell* cell = geometry::cells[cell_id].get();

          // Only consider cells which have a boundary condition.
          if (cell->vacuum_or_reflective() == false) continue;

          auto d_t = cell->distance_to_boundary_condition(pad.r_local, u_,
                                                          surface_token_);
          if (d_t.first < dist && std::abs(d_t.first - dist) > BOUNDRY_TOL) {
            double tmp_dist = d_t.first;
            int32_t tmp_token = std::abs(d_t.second);

            if (tmp_token) {
              dist = tmp_dist;
              token = tmp_token;
              surface_index = token - 1;
            } else {
              // Not actually a surface, so we continue
              continue;
            }

            btype = geometry::surfaces[static_cast<std::size_t>(surface_index)]
                        ->boundary();

            if (geometry::surfaces[static_cast<std::size_t>(surface_index)]
                    ->sign(pad.r_local, u_) < 0)
              token *= -1;
          }
        } else if (pad.type != GeoLilyPad::PadType::Cell) {
          auto uni_indx = universe_id_to_indx[pad.id];
          Universe* uni = geometry::universes[uni_indx].get();
          if (uni->has_boundary_conditions()) {
            Boundary uni_bound =
                uni->get_boundary_condition(pad.r_local, u_, surface_token_);
            if (uni_bound.distance < dist &&
                std::abs(uni_bound.distance - dist) > BOUNDRY_TOL) {
              dist = uni_bound.distance;
              token = uni_bound.token;
              surface_index = uni_bound.surface_index;
              btype = uni_bound.boundary_type;
            }
          }
        }
      }

      Boundary ret_bound(dist, surface_index, btype);
      ret_bound.token = token;

      // Distance to surface, and token, with sign indicating the positions
      // current orientation to the surface.
      return ret_bound;
    } else {
      // This else is mainly useful when plotting the geometry, and the window
      // extends beyond the defined geometry.
      return geometry::root_universe->get_boundary_condition(r_, u_,
                                                             surface_token_);
    }
  }

  Boundary get_nearest_boundary() const {
    if (!this->is_lost()) {
      // This must be first, so that we make sure that boundary condition
      // surfaces are given the higher priority. That is, we need to really be
      // closer to another surface, to ignore the boundary condition surface.
      auto bound = this->get_boundary_condition();

      double dist = bound.distance;
      BoundaryType btype = bound.boundary_type;
      int surface_index = bound.surface_index;
      int32_t token = bound.token;

      // Go up the entire tree
      for (const auto& pad : tree) {
        if (pad.type == GeoLilyPad::PadType::Lattice) {
          const auto lat_indx = universe_id_to_indx[pad.id];
          const auto& lat = geometry::universes[lat_indx];
          double d = lat->distance_to_tile_boundary(pad.r_local, u_, pad.tile);
          if (d < dist && std::abs(d - dist) > BOUNDRY_TOL) {
            dist = d;
            btype = BoundaryType::Normal;
            surface_index = -1;
            token = 0;
          }
        } else if (pad.type == GeoLilyPad::PadType::Cell) {
          auto cell_id = cell_id_to_indx[pad.id];
          auto d_t = geometry::cells[cell_id]->distance_to_boundary(
              pad.r_local, u_, surface_token_);
          if (d_t.first < dist && std::abs(d_t.first - dist) > BOUNDRY_TOL) {
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

      Boundary ret_bound(dist, surface_index, btype);
      ret_bound.token = token;

      // Distance to surface, and token, with sign indicating the positions
      // current orientation to the surface.
      return ret_bound;
    } else {
      return geometry::root_universe->lost_get_boundary(r_, u_, surface_token_);
    }
  }

  void cross_surface(Boundary d_t) {
    this->move(d_t.distance);
    // Use negative as we have changed orientation !
    surface_token_ = -d_t.token;
  }

  bool is_lost() const { return !current_cell; }

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
        fatal_error(mssg.str());
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
        const auto lat_indx = universe_id_to_indx[it->id];
        const auto& lat = geometry::universes[lat_indx];
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
      if (size == 0) return this->restart_get_current();
      tree.resize(size);

      // Now start at the last element, and get the new position.
      // We can get the info for the last universe, and then call get_cell from
      // that Universe to descend the geometry tree.
      auto uni_indx = universe_id_to_indx[tree.back().id];
      const auto& uni = geometry::universes[uni_indx];
      Position r_local = tree.back().r_local;
      tree.pop_back();
      current_cell = uni->get_cell(tree, r_local, u_, surface_token_);

      // If we couldn't get a cell, we need to call in the big guns, and
      // re-start from scratch.
      if (!current_cell) restart_get_current();

      // If we have a valid cell, we must now also fill the material.
      if (current_cell) {
        if (current_cell.cell->fill() == Cell::Fill::Universe) {
          fatal_error("Did not find a cell with a material.");
        }

        current_mat = current_cell.cell->material();
      }
    }
  }

  bool check_tree() const {
    return r_.x() == tree.front().r_local.x() &&
           r_.y() == tree.front().r_local.y() &&
           r_.z() == tree.front().r_local.z();
  }

  void do_reflection(Particle& p, Boundary boundary) {
    // Get the surface first
    if (boundary.surface_index < 0) {
      fatal_error("Bad surface index in Tracker::do_reflection");
    }
    const std::shared_ptr<Surface>& surface =
        geometry::surfaces[static_cast<std::size_t>(boundary.surface_index)];

    int32_t token = geometry::id_to_token(static_cast<int32_t>(surface->id()));
    if (surface->sign(p.r(), p.u()) < 0) token *= -1;
    this->set_surface_token(token);

    // Get new Position object to temporarily contain the current position
    Position r_pre_refs = p.r();

    // In the event a reflection occurs immediately after another, we must
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
  UniqueCell current_cell;
  // Token for current surface which particle is on. The index for the
  // surface in the surface vector is std::abs(token) - 1. Token is NOT
  // the surface id !!
  int32_t surface_token_ = 0;

};  // Tracker

#endif  // MG_TRACKER_H
