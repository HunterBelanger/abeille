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
#ifndef APPROXIMATE_MESH_CANCELATOR_H
#define APPROXIMATE_MESH_CANCELATOR_H

#include <simulation/cancelator.hpp>
#include <unordered_map>

class ApproximateMeshCancelator : public Cancelator {
 public:
  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz);
  ApproximateMeshCancelator(Position low, Position hi, uint32_t Nx, uint32_t Ny,
                            uint32_t Nz, std::vector<double> energy_bounds);

  bool add_particle(BankedParticle& p) override final;
  void perform_cancellation(pcg32& rng) override final;
  std::vector<BankedParticle> get_new_particles(pcg32& rng) override final;
  void clear() override final;

 private:
  Position r_low, r_hi;
  std::vector<double> energy_edges;
  std::array<uint32_t, 4> shape;
  double dx, dy, dz;
  std::unordered_map<int, std::vector<BankedParticle*>> bins;
};

std::shared_ptr<ApproximateMeshCancelator> make_approximate_mesh_cancelator(
    const YAML::Node& node);

#endif
