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
