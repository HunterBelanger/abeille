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
#ifndef GEO_LILY_PAD_H
#define GEO_LILY_PAD_H

#include <utils/position.hpp>

#include <array>
#include <cstdint>

struct GeoLilyPad {
  enum class PadType { Universe, Lattice, Cell };

  GeoLilyPad() {}
  GeoLilyPad(PadType t, uint32_t i, Position r, std::array<int32_t, 3> ti,
             bool ou)
      : type(t), id(i), r_local(r), tile(ti), in_lattice_outside_universe(ou) {}

  PadType type = PadType::Universe;
  uint32_t id = 0;

  // r_local is the position within the lattice,
  // NOT the position within the tile !!
  Position r_local = Position();
  std::array<int32_t, 3> tile = {0, 0, 0};
  bool in_lattice_outside_universe = false;
};

#endif
