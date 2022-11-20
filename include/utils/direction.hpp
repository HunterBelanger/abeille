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
