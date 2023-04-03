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
#include <simulation/approximate_mesh_cancelator.hpp>
#include <simulation/basic_exact_mg_cancelator.hpp>
#include <simulation/cancelator.hpp>
#include <simulation/exact_mg_cancelator.hpp>
#include <utils/error.hpp>
#include <utils/settings.hpp>

std::shared_ptr<Cancelator> make_cancelator(const YAML::Node& node) {
  if (!node["type"] || !node["type"].IsScalar()) {
    fatal_error("Invalid type entry for cancelator.");
  }
  std::string type = node["type"].as<std::string>();

  std::shared_ptr<Cancelator> cancelator = nullptr;

  if (type == "approximate") {
    cancelator = make_approximate_mesh_cancelator(node);
  } else if (type == "basic-exact") {
    // Check that we are not using surface-tracking !
    if (settings::tracking == settings::TrackingMode::SURFACE_TRACKING) {
      fatal_error(
          "basic-exect cancelators may not be used with surface-tracking.");
    }

    // Check that we are not in CE mode !
    if (settings::energy_mode == settings::EnergyMode::CE) {
      fatal_error(
          "basic-exact cancelators are not currently supported for continuous "
          "energy.");
    }

    // looks like we meet all requirments, make the cancelator
    cancelator = make_basic_exact_mg_cancelator(node);
  } else if (type == "exact") {
    // Check that we are not using surface-tracking !
    if (settings::tracking == settings::TrackingMode::SURFACE_TRACKING) {
      fatal_error("exect cancelators may not be used with surface-tracking.");
    }

    // Check that we are not in CE mode !
    if (settings::energy_mode == settings::EnergyMode::CE) {
      fatal_error(
          "exact cancelators are not currently supported for continuous "
          "energy.");
    }

    // looks like we meet all requirments, make the cancelator
    cancelator = make_exact_mg_cancelator(node);
  } else {
    fatal_error("Unkown cancelator type \"" + type + "\".");
  }

  return cancelator;
}
