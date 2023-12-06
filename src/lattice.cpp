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
#include <geometry/geometry.hpp>
#include <geometry/lattice.hpp>

#include <cstdint>

Lattice::Lattice(uint32_t i_id, std::string i_name)
    : lattice_universes{}, outer_universe_index{-1}, id_{i_id}, name_{i_name} {}

void Lattice::set_outisde_universe(int32_t univ) {
  outer_universe_index = univ;
}

Universe* Lattice::outer_universe() const {
  if (outer_universe_index < 0) return nullptr;

  return geometry::universes[static_cast<std::size_t>(outer_universe_index)]
      .get();
}

const Universe* Lattice::get_universe(std::size_t ind) const {
  if (lattice_universes[ind] < 0) return nullptr;

  std::size_t uni_indx = static_cast<std::size_t>(lattice_universes[ind]);

  return geometry::universes[uni_indx].get();
}

 Universe* Lattice::get_universe(std::size_t ind)  {
  if (lattice_universes[ind] < 0) return nullptr;

  std::size_t uni_indx = static_cast<std::size_t>(lattice_universes[ind]);

  return geometry::universes[uni_indx].get();
}

uint32_t Lattice::id() const { return id_; }

std::string Lattice::name() const { return name_; }
