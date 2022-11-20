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
#ifndef TRACK_LENGTH_MESH_TALLY_H
#define TRACK_LENGTH_MESH_TALLY_H

#include <yaml-cpp/yaml.h>

#include <simulation/mesh_tally.hpp>

class TrackLengthMeshTally : public MeshTally {
 public:
  TrackLengthMeshTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                       uint64_t nz, const std::vector<double>& ebounds,
                       Quantity q, std::string fname, uint32_t mt = 0)
      : MeshTally(low, hi, nx, ny, nz, ebounds, fname), quantity(q), mt(mt) {}

  void score_flight(const Particle& p, double d, MaterialHelper& mat);

 private:
  Quantity quantity;
  uint32_t mt;

  void initialize_indices(const Position& r, const Direction& u, int& i, int& j,
                          int& k, std::array<int, 3>& on);
  std::pair<double, int> distance_to_next_index(const Position& r,
                                                const Direction& u, int i,
                                                int j, int k,
                                                const std::array<int, 3>& on);
  bool find_entry_point(Position& r, const Direction& u,
                        double& d_flight) const;
  void update_indices(int key, int& i, int& j, int& k, std::array<int, 3>& on);

  double get_base_score(const Particle& p, MaterialHelper& mat) const;

  std::string quantity_str() const override final;
};

std::shared_ptr<TrackLengthMeshTally> make_track_length_mesh_tally(
    const YAML::Node& node);

#endif
