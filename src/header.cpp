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
#include <utils/constants.hpp>
#include <utils/header.hpp>

const std::string logo =
    "\n"
    "                    ______   __                  __  __  __           \n"
    "                   /      \\ |  \\                |  \\|  \\|  \\         "
    " \n"
    "                  |  $$$$$$\\| $$____    ______   \\$$| $$| $$  ______  \n"
    "                  | $$__| $$| $$    \\  /      \\ |  \\| $$| $$ /      \\ "
    "\n"
    "                  | $$    $$| $$$$$$$\\|  $$$$$$\\| $$| $$| $$|  $$$$$$\n"
    "                  | $$$$$$$$| $$  | $$| $$    $$| $$| $$| $$| $$    $$\n"
    "                  | $$  | $$| $$__/ $$| $$$$$$$$| $$| $$| $$| $$$$$$$$\n"
    "                  | $$  | $$| $$    $$ \\$$     \\| $$| $$| $$ \\$$     \n"
    "                   \\$$   \\$$ \\$$$$$$$   \\$$$$$$$ \\$$ \\$$ \\$$  "
    "\\$$$$$$$\n"
    "\n\n"
    "                          .-=-=-=-.        \n"
    "                        (`-=-=-=-=-`)      \n"
    "                      (`-=-=-=-=-=-=-`)                      __\n"
    "                     (`-=-=-=-=-=-=-=-`)                    // \\  \n"
    "                    (`-=-=-=-( @ )-=-=-`)                   \\\\_/ //  \n"
    "                    (`-=-=-=-=-=-=-=-=-`) ''-.._.-''-.._.. -(||)(') \n"
    "                    (`-=-=-=-=-=-=-=-=-`)                   ''' \n"
    "                    (`-=-=-=-=-=-=-=-=-`)  \n"
    "                    (`-=-=-=-=-=-=-=-=-`)  \n"
    "                     (`-=-=-=-=-=-=-=-`)   \n"
    "                      (`-=-=-=-=-=-=-`)    \n"
    "                        (`-=-=-=-=-`)      \n"
    "                         `-=-=-=-=-`       \n"
    "\n\n";

const std::string header =
    "                              A Monte Carlo Transport Code\n\n";

const std::string info =
    " Copyright (C) 2019-2023, Hunter Belanger\n"
    " Copyright (C) 2021-2022, Commissariat à l'énergie atomique et aux "
    "énergies alternatives (CEA)\n"
    " Released under the terms and conditions of the GPLv3 license\n"
    " Written by Hunter Belanger\n";

const std::string help =
    " Usage:\n"
#ifdef ABEILLE_USE_OMP
    "   abeille (--input FILE) [--threads NUM --output FILE]\n"
#else
    "   abeille (--input FILE) [--output FILE]\n"
#endif
    "   abeille (--input FILE) (--gui-plot | --plot) [--threads NUM]\n"
    "   abeille (-h | --help)\n"
    "   abeille (-v | --version)\n\n"

    " Options:\n"
    "   -h --help         Show this help message\n"
    "   -v --version      Show version number\n"
    "   -i --input FILE   Set input file\n"
#ifdef ABEILLE_USE_OMP
    "   -t --threads NUM  Set number of OpenMP threads\n"
#endif
    "   -o --output FILE  Set output file\n"
#ifdef ABEILLE_GUI_PLOT
    "   -g --gui-plot     Interactive GUI plotter\n"
#endif
    "   -p --plot         Generate plots from input file\n";

const std::string version_string =
    " Abeille " ABEILLE_VERSION_STRING
    " "
#if defined(DEVELOPMENT_VERSION)
    "(Development)"
#endif
    "\n"
    " Git Hash " ABEILLE_GIT_HASH
    "\n"
    " Compiler " ABEILLE_COMPILER_NAME " " ABEILLE_COMPILER_VERSION
    "\n"
#if defined(ABEILLE_GUI_PLOT)
    " Compiled with GUI plotter.\n"
#endif
#if defined(ABEILLE_USE_OMP)
    " Compiled with OpenMP.\n"
#endif
#if defined(ABEILLE_USE_MPI)
    " Compiled with MPI.\n"
#endif
    "";
