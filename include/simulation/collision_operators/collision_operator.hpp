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
#ifndef COLLISION_OPERATOR_H
#define COLLISION_OPERATOR_H

#include <materials/material_helper.hpp>
#include <simulation/particle.hpp>
#include <tallies/tallies.hpp>

#include <concepts>
#include <string>

template <typename C>
concept CollisionOperator = requires(C c, Particle& p, MaterialHelper& mat,
                                     ThreadLocalScores& thread_scores,
                                     const std::string& base) {
  { c.collision(p, mat, thread_scores) } -> std::same_as<void>;
  { c.write_output_info(base) } -> std::same_as<void>;
};

#endif