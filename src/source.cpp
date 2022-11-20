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
#include <simulation/source.hpp>
#include <simulation/tracker.hpp>
#include <utils/error.hpp>

Source::Source(std::shared_ptr<SpatialDistribution> spatial,
               std::shared_ptr<DirectionDistribution> direction,
               std::shared_ptr<EnergyDistribution> energy, bool fissile_only,
               double weight)
    : spatial_(spatial),
      direction_(direction),
      energy_(energy),
      fissile_only_(fissile_only),
      weight_(weight) {
  if (weight_ <= 0.) {
    std::string mssg = "Source weight must be greater than zero.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

Particle Source::generate_particle(pcg32& rng) const {
  // Sample Direction
  Direction u = direction_->sample(rng);

  // Sample energy
  double E = energy_->sample(rng);

  // Sample position
  Position r = spatial_->sample(rng);
  Tracker trkr(r, u);

  // Keep sampling new position untill a position inside the geometry
  // is found.
  while (trkr.is_lost()) {
    r = spatial_->sample(rng);
    trkr.set_r(r);
    trkr.restart_get_current();
  }

  // If we only want positions inside of fissile regions, then we need
  // to check this.
  if (fissile_only_) {
    int count = 0;
    while (trkr.material()->fissile() == false || trkr.is_lost()) {
      if (count == 201) {
        std::string mssg =
            "Exceded 200 samplings of position in fissile-only source.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      r = spatial_->sample(rng);
      trkr.set_r(r);
      trkr.restart_get_current();
      count++;
    }
  }

  return Particle(r, u, E, 1.0);
}

std::shared_ptr<Source> make_source(YAML::Node source_node) {
  // Make sure it is a map
  if (!source_node.IsMap()) {
    std::string mssg = "Source entry must be a map.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get weight
  if (!source_node["weight"] || !source_node["weight"].IsScalar()) {
    std::string mssg = "No weight given to source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  double wgt = source_node["weight"].as<double>();

  // Get fissile only
  bool fissile_only = false;
  if (source_node["fissile-only"] && source_node["fissile-only"].IsScalar()) {
    fissile_only = source_node["fissile-only"].as<bool>();
  } else if (source_node["fissile-only"]) {
    std::string mssg = "Invalid fissile-only entry in source description.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get spatial distribution
  if (!source_node["spatial"] || !source_node["spatial"].IsMap()) {
    std::string mssg =
        "No valid spatial distribution entry provided for source distribution.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::shared_ptr<SpatialDistribution> spatial =
      make_spatial_distribution(source_node["spatial"]);

  // Get direction distribution
  if (!source_node["direction"] || !source_node["direction"].IsMap()) {
    std::string mssg =
        "No valid direction distribution entry provided for source "
        "distribution.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::shared_ptr<DirectionDistribution> direction =
      make_direction_distribution(source_node["direction"]);

  // Get energy distribution
  if (!source_node["energy"] || !source_node["energy"].IsMap()) {
    std::string mssg =
        "No valid energy distribution entry provided for source distribution.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::shared_ptr<EnergyDistribution> energy =
      make_energy_distribution(source_node["energy"]);

  return std::make_shared<Source>(spatial, direction, energy, fissile_only,
                                  wgt);
}
