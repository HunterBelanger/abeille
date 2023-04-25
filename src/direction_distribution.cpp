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
#include <memory>
#include <simulation/cone.hpp>
#include <simulation/direction_distribution.hpp>
#include <simulation/isotropic.hpp>
#include <simulation/mono_directional.hpp>
#include <utils/error.hpp>

std::shared_ptr<DirectionDistribution> make_direction_distribution(
    YAML::Node node) {
  // Get the type of the distribution
  if (!node["type"] || !node["type"].IsScalar()) {
    fatal_error("No valid type provided to direction distribution entry.");
  }

  std::string type = node["type"].as<std::string>();

  std::shared_ptr<DirectionDistribution> dist = nullptr;

  if (type == "mono-directional") {
    dist = make_mono_directional_distribution(node);
  } else if (type == "isotropic") {
    dist = std::make_shared<Isotropic>();
  } else if (type == "cone") {
    dist = make_cone_distribution(node);
  } else {
    fatal_error("Invalid direction distribution type " + type + ".");
  }

  return dist;
}
