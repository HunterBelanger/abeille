/*
 * Abeille Monte Carlo Code
 * Copyright 2019-2023, Hunter Belanger
 * Copyright 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Abeille Monte Carlo code (Abeille).
 *
 * Abeille is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Abeille is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Abeille. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

//============================================================================
// Versioning
constexpr int VERSION_MAJOR{0};
constexpr int VERSION_MINOR{4};
constexpr int VERSION_PATCH{0};
constexpr int COPYRIGHT_YEAR{2023};
#define DEVELOPMENT_VERSION
#define ABEILLE_VERSION_STRING "0.4.0"

//============================================================================
// Mathematical and Physical Constants
constexpr double INF{std::numeric_limits<double>::max()};
constexpr int INF_INT{std::numeric_limits<int>::max()};
constexpr double PI{3.14159265358979323846264338327950288};
constexpr double EV_TO_K{1.160451812E4};
constexpr double MEV_TO_EV{1.E6};
constexpr double EV_TO_MEV{1.E-6};
constexpr double N_MASS_EV{939.56542052 * MEV_TO_EV};
constexpr double N_MASS_AMU{1.00866491595};
constexpr double C_CM_S{29979245800.0};
constexpr double N_AVAGADRO{0.6022140857};  // [10^24 / mol]
constexpr double TOLERANCE{0.0001};

//============================================================================
// Program Parameters
constexpr double SURFACE_COINCIDENT{1E-12};
constexpr double BOUNDRY_TOL{500. * SURFACE_COINCIDENT};

#endif
