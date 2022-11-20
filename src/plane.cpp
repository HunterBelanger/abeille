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
#include <geometry/surfaces/plane.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

Plane::Plane(double A_, double B_, double C_, double D_, BoundaryType bound,
             uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, A{A_}, B{B_}, C{C_}, D{D_} {}

int Plane::sign(const Position& r, const Direction& u) const {
  double eval = A * r.x() + B * r.y() + C * r.z() - D;
  if (eval > SURFACE_COINCIDENT)
    return 1;
  else if (eval < -SURFACE_COINCIDENT)
    return -1;
  else {
    if (u.dot(norm(r)) > 0.)
      return 1;
    else
      return -1;
  }
}

double Plane::distance(const Position& r, const Direction& u,
                       bool on_surf) const {
  double num = D - A * r.x() - B * r.y() - C * r.z();
  double denom = A * u.x() + B * u.y() + C * u.z();
  double d = num / denom;
  if (on_surf || std::abs(d) < SURFACE_COINCIDENT || denom == 0.)
    return INF;
  else if (d < 0.)
    return INF;
  else
    return d;
}

Direction Plane::norm(const Position& /*r*/) const { return {A, B, C}; }

//===========================================================================
// Non-Member Functions
std::shared_ptr<Plane> make_plane(YAML::Node surface_node) {
  // Variables for surface
  double A = 0., B = 0., C = 0., D = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get A
  if (surface_node["A"])
    A = surface_node["A"].as<double>();
  else {
    std::string mssg = "Plane surface must have A defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get B
  if (surface_node["B"])
    B = surface_node["B"].as<double>();
  else {
    std::string mssg = "Plane surface must have B defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get C
  if (surface_node["C"])
    C = surface_node["C"].as<double>();
  else {
    std::string mssg = "Plane surface must have C defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get D
  if (surface_node["D"])
    D = surface_node["D"].as<double>();
  else {
    std::string mssg = "Plane surface must have D defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get boundary type
  if (surface_node["boundary"]) {
    std::string boundary_string = surface_node["boundary"].as<std::string>();
    if (boundary_string == "vacuum")
      boundary = BoundaryType::Vacuum;
    else if (boundary_string == "reflective")
      boundary = BoundaryType::Reflective;
    else if (boundary_string == "normal")
      boundary = BoundaryType::Normal;
    else {
      std::string mssg = "Unknown boundary type \"" + boundary_string + "\".";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    boundary = BoundaryType::Normal;
  }

  // Get id
  if (surface_node["id"])
    id = surface_node["id"].as<uint32_t>();
  else {
    std::string mssg =
        "Surface must have an id attribute with a unique"
        " positive integer.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get name
  if (surface_node["name"])
    name = surface_node["name"].as<std::string>();
  else
    name = "";

  // Make and return surface
  return std::make_shared<Plane>(A, B, C, D, boundary, id, name);
}
