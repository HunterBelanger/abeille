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
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <geometry/boundary.hpp>
#include <geometry/cell.hpp>
#include <geometry/geo_lily_pad.hpp>
#include <geometry/lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <geometry/universe.hpp>

#include <cstddef>

namespace geometry {

//==========================================================================
// Vectors for Geometry Objects
//--------------------------------------------------------------------------
// All surfaces in problem
extern std::vector<std::shared_ptr<Surface>> surfaces;

// All cells in problem
extern std::vector<std::shared_ptr<Cell>> cells;

// Root Universe containing problem geometry, may be a
// CellUniverse or a LatticeUniverse
extern std::shared_ptr<Universe> root_universe;

// All universes in problem(including root universe at 0)
extern std::vector<std::shared_ptr<Universe>> universes;

//==========================================================================
// Functions
UniqueCell get_cell(const Position& r, const Direction& u, int32_t on_surf = 0);

UniqueCell get_cell(std::vector<GeoLilyPad>& stack, const Position& r,
                    const Direction& u, int32_t on_surf = 0);

int32_t id_to_token(int32_t id);

void do_reflection(Particle& p, Boundary boundary);

}  // namespace geometry

#endif
