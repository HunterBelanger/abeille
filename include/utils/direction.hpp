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
#ifndef DIRECTION_H
#define DIRECTION_H

#include <utils/constants.hpp>
#include <utils/vector.hpp>

//============================================================================
// Direction Class
//----------------------------------------------------------------------------
class Direction : public Vector {
 public:
  Direction() : Vector(1., 0., 0.) {}
  Direction(double i_x, double i_y, double i_z) : Vector{i_x, i_y, i_z} {
    double magnitude = this->norm();
    this->x_ /= magnitude;
    this->y_ /= magnitude;
    this->z_ /= magnitude;
  }
  Direction(double mu, double phi) : Vector(0., 0., 1.) {
    if (mu < -1.)
      mu = -1.;
    else if (mu > 1.)
      mu = 1.;

    if (phi < 0.)
      phi = 0.;
    else if (phi > 2 * PI)
      phi = 2 * PI;

    this->x_ = std::sqrt(1. - mu * mu) * std::cos(phi);
    this->y_ = std::sqrt(1. - mu * mu) * std::sin(phi);
    this->z_ = mu;

    double magnitude = this->norm();
    this->x_ /= magnitude;
    this->y_ /= magnitude;
    this->z_ /= magnitude;
  }
  ~Direction() = default;

};  // Direction

//============================================================================
// Addition Operators
inline Vector operator+(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() + d2.x(), d1.y() + d2.y(), d1.z() + d2.z());
}

inline Vector operator+(const Direction& d, const Vector& v) {
  return Vector(d.x() + v.x(), d.y() + v.y(), d.z() + v.z());
}

inline Vector operator+(const Vector& v, const Direction& d) {
  return Vector(d.x() + v.x(), d.y() + v.y(), d.z() + v.z());
}

//============================================================================
// Subtraction Operators
inline Vector operator-(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() - d2.x(), d1.y() - d2.y(), d1.z() - d2.z());
}

inline Vector operator-(const Direction& d, const Vector& v) {
  return Vector(d.x() - v.x(), d.y() - v.y(), d.z() - v.z());
}

inline Vector operator-(const Vector& v, const Direction& d) {
  return Vector(v.x() - d.x(), v.y() - d.y(), v.z() - d.z());
}

//============================================================================
// Dot Product Operators
inline double operator*(const Direction& d, const Direction& v) {
  return d.dot(v);
}

inline double operator*(const Direction& d, const Vector& v) {
  return d.dot(v);
}

inline double operator*(const Vector& v, const Direction& d) {
  return d.dot(v);
}

//============================================================================
// Scaling Operators
inline Vector operator*(const Direction& d, double c) {
  return Vector(d.x() * c, d.y() * c, d.z() * c);
}

inline Vector operator*(double c, const Direction& d) {
  return Vector(d.x() * c, d.y() * c, d.z() * c);
}

inline Vector operator/(const Direction& d, double c) {
  return Vector(d.x() / c, d.y() / c, d.z() / c);
}

inline std::ostream& operator<<(std::ostream& output, const Direction& d) {
  output << "<<" << d.x() << "," << d.y() << "," << d.z() << ">>";
  return output;
}

// This function returneds the direction u, rotated about the
// axis of rotation defined by a x b.
inline Direction rotate_direction(Direction u, double mu, double phi) {
  double cos = std::cos(phi);
  double sin = std::sin(phi);
  double sqrt_mu = std::sqrt(1. - mu * mu);
  double sqrt_w = std::sqrt(1. - u.z() * u.z());

  double ux, uy, uz;

  if (sqrt_w > 1.E-10) {
    ux = mu * u.x() + sqrt_mu * (u.x() * u.z() * cos - u.y() * sin) / sqrt_w;
    uy = mu * u.y() + sqrt_mu * (u.y() * u.z() * cos + u.x() * sin) / sqrt_w;
    uz = mu * u.z() - sqrt_mu * sqrt_w * cos;
  } else {
    double sqrt_v = std::sqrt(1. - u.y() * u.y());

    ux = mu * u.x() + sqrt_mu * (u.x() * u.y() * cos + u.z() * sin) / sqrt_v;
    uy = mu * u.y() - sqrt_mu * sqrt_v * cos;
    uz = mu * u.z() + sqrt_mu * (u.y() * u.z() * cos - u.x() * sin) / sqrt_v;
  }

  return {ux, uy, uz};
}

#endif
