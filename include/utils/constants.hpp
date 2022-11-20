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
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <limits>

//============================================================================
// Versioning
constexpr int VERSION_MAJOR{0};
constexpr int VERSION_MINOR{3};
constexpr int VERSION_PATCH{0};
constexpr bool DEVELOPMENT_VERSION{false};
//#define DEV_VERSION  // TO BE REMOVED FOR NON-DEV VERSIONS
constexpr int COPYRIGHT_YEAR{2022};
#define ABEILLE_VERSION_STRING "0.3.0"

//============================================================================
// Mathematical and Physical Constants
constexpr double INF{std::numeric_limits<double>::max()};
constexpr int INF_INT{std::numeric_limits<int>::max()};
constexpr double PI{3.14159265358979323846264338327950288};
constexpr double EV_TO_K{1.160451812E4};
constexpr double MEV_TO_EV{1.E6};
constexpr double EV_TO_MEV{1.E-6};
constexpr double N_MASS_EV{939.56542052 * MEV_TO_EV};
constexpr double C_CM_S{29979245800.0};
constexpr double N_AVAGADRO{0.6022140857};  // [10^24 / mol]
constexpr double TOLERANCE{0.0001};

//============================================================================
// Program Parameters
constexpr double SURFACE_COINCIDENT{1E-12};

#endif
