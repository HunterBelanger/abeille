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
#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

//============================================================================
// Vector Class
//----------------------------------------------------------------------------
class Vector {
 public:
  Vector(double i_x, double i_y, double i_z) : x_{i_x}, y_{i_y}, z_{i_z} {};
  ~Vector() = default;

  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }

  double dot(const Vector& v) const {
    return x_ * v.x() + y_ * v.y() + z_ * v.z();
  }

  double norm() const { return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_); }

 protected:
  double x_;
  double y_;
  double z_;

};  // Vector

//============================================================================
// Overloaded Operator Declarations
//----------------------------------------------------------------------------
inline Vector operator+(const Vector& v1, const Vector& v2) {
  return Vector(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline Vector operator-(const Vector& v1, const Vector& v2) {
  return Vector(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

inline Vector operator*(const Vector& v, double d) {
  return Vector(v.x() * d, v.y() * d, v.z() * d);
}

inline Vector operator*(double d, const Vector& v) { return v * d; }

inline Vector operator/(const Vector& v, double d) {
  return Vector(v.x() / d, v.y() / d, v.z() / d);
}

inline double operator*(const Vector& v1, const Vector& v2) {
  return v1.dot(v2);
}

inline std::ostream& operator<<(std::ostream& output, const Vector& v) {
  output << "<" << v.x() << "," << v.y() << "," << v.z() << ">";
  return output;
}

#endif
