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
#ifndef POSITION_H
#define POSITION_H

#include <utils/vector.hpp>

//============================================================================
// Position Class
//----------------------------------------------------------------------------
class Position : public Vector {
 public:
  Position() : Vector{0., 0., 0.} {}
  Position(double i_x, double i_y, double i_z) : Vector{i_x, i_y, i_z} {};
  ~Position() = default;
};  // Position

//============================================================================
// Comparison Operators
inline bool operator==(const Position& p1, const Position& p2) {
  return (p1.x() == p2.x()) && (p1.y() == p2.y()) && (p1.z() == p2.z());
}

inline bool operator!=(const Position& p1, const Position& p2) {
  return !(p1 == p2);
}

//============================================================================
// Addition Operators
inline Position operator+(const Position& p1, const Position& p2) {
  return Position(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}

inline Position operator+(const Vector& v1, const Position& p2) {
  return Position(v1.x() + p2.x(), v1.y() + p2.y(), v1.z() + p2.z());
}

inline Position operator+(const Position& p2, const Vector& v1) {
  return Position(v1.x() + p2.x(), v1.y() + p2.y(), v1.z() + p2.z());
}

//============================================================================
// Subtraction Operators
inline Position operator-(const Position& p1, const Position& p2) {
  return Position(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}

inline Position operator-(const Vector& v1, const Position& p2) {
  return Position(v1.x() - p2.x(), v1.y() - p2.y(), v1.z() - p2.z());
}

inline Position operator-(const Position& p2, const Vector& v1) {
  return Position(p2.x() - v1.x(), p2.y() - v1.y(), p2.z() - v1.z());
}

//============================================================================
// Dot Products
inline double operator*(const Position& p1, const Position& p2) {
  return p1.dot(p2);
}

inline double operator*(const Vector& p1, const Position& p2) {
  return p1.dot(p2);
}

inline double operator*(const Position& p2, const Vector& p1) {
  return p1.dot(p2);
}

inline Position operator*(const Position& p, double d) {
  return Position(p.x() * d, p.y() * d, p.z() * d);
}

inline Position operator*(double d, const Position& p) { return p * d; }

inline Position operator/(const Position& p, double d) {
  return Position(p.x() / d, p.y() / d, p.z() / d);
}

inline std::ostream& operator<<(std::ostream& output, const Position& p) {
  output << "(" << p.x() << "," << p.y() << "," << p.z() << ")";
  return output;
}

#endif
