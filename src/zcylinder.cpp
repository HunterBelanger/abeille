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
#include <geometry/surfaces/zcylinder.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

ZCylinder::ZCylinder(double x_, double y_, double r_, BoundaryType bound,
                     uint32_t i_id, std::string i_name)
    : Surface{bound, i_id, i_name}, x0{x_}, y0{y_}, R{r_} {}

int ZCylinder::sign(const Position& r, const Direction& u) const {
  double x = r.x() - x0;
  double y = r.y() - y0;
  double eval = y * y + x * x - R * R;
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

double ZCylinder::distance(const Position& r, const Direction& u,
                           bool on_surf) const {
  double a = u.y() * u.y() + u.x() * u.x();
  if (a == 0.) return INF;

  double x = r.x() - x0;
  double y = r.y() - y0;
  double k = y * u.y() + x * u.x();
  double c = y * y + x * x - R * R;
  double quad = k * k - a * c;

  if (quad < 0.)
    return INF;
  else if (on_surf || std::abs(c) < SURFACE_COINCIDENT) {
    if (k >= 0.)
      return INF;
    else
      return (-k + std::sqrt(quad)) / a;
  } else if (c < 0.) {
    return (-k + std::sqrt(quad)) / a;
  } else {
    double d = (-k - std::sqrt(quad)) / a;
    if (d < 0.)
      return INF;
    else
      return d;
  }
}

Direction ZCylinder::norm(const Position& r) const {
  double x = r.x() - x0;
  double y = r.y() - y0;
  return {x, y, 0.};
}

//===========================================================================
// Non-Member Functions
std::shared_ptr<ZCylinder> make_zcylinder(YAML::Node surface_node) {
  // Variables for surface
  double x0 = 0., y0 = 0., r = 0.;
  BoundaryType boundary = BoundaryType::Normal;
  uint32_t id = 1;
  std::string name = "";

  // Get x0
  if (surface_node["x0"])
    x0 = surface_node["x0"].as<double>();
  else {
    std::string mssg = "ZCylinder surface must have x0 defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get y0
  if (surface_node["y0"])
    y0 = surface_node["y0"].as<double>();
  else {
    std::string mssg = "ZCylinder surface must have y0 defined.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get r
  if (surface_node["r"])
    r = surface_node["r"].as<double>();
  else {
    std::string mssg = "ZCylinder surface must have r defined.";
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
  return std::make_shared<ZCylinder>(x0, y0, r, boundary, id, name);
}
