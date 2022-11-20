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
