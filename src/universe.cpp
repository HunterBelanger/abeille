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
#include <geometry/universe.hpp>
#include <utils/constants.hpp>

Universe::Universe(uint32_t i_id, std::string i_name)
    : cell_offset_map{}, id_{i_id}, name_{i_name} {}

std::array<int32_t, 3> Universe::get_tile(Position /*r*/,
                                          Direction /*u*/) const {
  return {0, 0, 0};
}

double Universe::distance_to_tile_boundary(
    Position /*r_local*/, Direction /*u*/,
    std::array<int32_t, 3> /*tile*/) const {
  return INF;
}

uint32_t Universe::id() const { return id_; }

std::string Universe::name() const { return name_; }
