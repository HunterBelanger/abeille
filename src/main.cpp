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
#include <docopt.h>

#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <plotting/plotter.hpp>
#include <stdexcept>
#include <string>
#include <utils/constants.hpp>
#include <utils/header.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>
#include <utils/timer.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

void signal_handler(int signal) {
  if (!simulation->signaled) {
    // Signal has been recieved. Set the signaled
    // flag in simulation. It will stop after the
    // current generation's particles have all been
    // transported as we can't break inside an
    // OpenMP loop.
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      simulation->signaled = true;
      if (signal == SIGTERM) simulation->terminate = true;
    }
  } else {
    // A second signal has been recieved. We shall obey the master's
    // wishes, and blow up everything, stoping the simulation imediately.
    // No tallies will be saved, and the program with exit with an error
    // code of 1.
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      Output::instance()->write("\n !! >>> FORCE KILLED BY USER <<< !!\n");
      mpi::abort_mpi();
      std::exit(1);
    }
  }
}

bool exists(std::string fname) {
  std::ifstream file(fname);
  return file.good();
}

int main(int argc, char** argv) {
  settings::alpha_omega_timer.start();

  mpi::initialize_mpi(&argc, &argv);

  std::atexit(mpi::finalize_mpi);

  // Make help message string for docopt
  std::string help_message = version_string + "\n" + info + "\n" + help;

  // Initialize docopt
  std::map<std::string, docopt::value> args =
      docopt::docopt(help_message, {argv + 1, argv + argc}, true,
                     version_string + "\n" + info);

  std::string output_filename;
  if (args["--output"])
    output_filename = args["--output"].asString();
  else
    output_filename = "output.txt";
  Output::set_output_filename(output_filename);
  mpi::synchronize();

#ifdef ABEILLE_GUI_PLOT
  if (args["--gui-plot"].asBool() || args["--plot"].asBool()) {
#else
  if (args["--plot"].asBool()) {
#endif
    // If using OpenMP, get number of threads requested
#ifdef _OPENMP
    int num_omp_threads;
    if (args["--threads"])
      num_omp_threads = std::stoi(args["--threads"].asString());
    else
      num_omp_threads = omp_get_max_threads();

    // If number of threads requested greater than system max, use system max
    if (omp_get_max_threads() < num_omp_threads) {
      num_omp_threads = omp_get_max_threads();
    } else {
      omp_set_num_threads(num_omp_threads);
    }
#endif
    // Get input file
    std::string input_fname = "";
    if (exists({args["--input"].asString()}))
      input_fname = args["--input"].asString();
    else if (exists("input.yaml"))
      input_fname = "input.yaml";

    if (input_fname.size() > 0) {
      print_header();
      Output::instance()->write("\n");

      try {
        // Begin plotting system
#ifdef ABEILLE_GUI_PLOT
        if (args["--gui-plot"].asBool()) {
          plotter::gui(input_fname);
        } else {
          plotter::plotter(input_fname);
        }
#else
        plotter::plotter(input_fname);
#endif
      } catch (const std::runtime_error& err) {
        std::string mssg = err.what();
        Output::instance()->write(" FATAL ERROR: " + mssg + ".\n");
        std::exit(1);
      }

      // After done with plotter, exit program
      return 0;
    } else {
      std::cout << " ERROR: Input file not found.\n";
    }
  } else if (exists({args["--input"].asString()})) {  // Run program
    std::string input_filename;

    if (args["--input"])
      input_filename = args["--input"].asString();
    else
      input_filename = "input.txt";

      // If using OpenMP, get number of threads requested
#ifdef _OPENMP
    int num_omp_threads;
    if (args["--threads"])
      num_omp_threads = std::stoi(args["--threads"].asString());
    else
      num_omp_threads = omp_get_max_threads();

    // If number of threads requested greater than system max, use system max
    if (omp_get_max_threads() < num_omp_threads) {
      num_omp_threads = omp_get_max_threads();
    } else {
      omp_set_num_threads(num_omp_threads);
    }
#endif

    // Print output header if rank 0
    mpi::synchronize();
    print_header();
    Output::instance()->write("\n");

    // Parse input file
    bool parsed_file = false;
    try {
      parse_input_file(args["--input"].asString());
      parsed_file = true;
    } catch (const std::runtime_error& err) {
      std::string mssg = err.what();
      Output::instance()->write(" FATAL ERROR: " + mssg + ".\n");
      parsed_file = false;
    }

    if (parsed_file) {
      // Try to run simulation
      simulation->initialize();

      // Setup signal catching
      std::signal(SIGINT, signal_handler);
      std::signal(SIGTERM, signal_handler);

      simulation->run();
    }

  } else {
    std::cout << " ERROR: Input file not found.\n";
  }

#ifdef ABEILLE_USE_MPI
  if (mpi::size > 1) {
    Output::instance()->write(
        " Total MPI time: " + std::to_string(mpi::timer.elapsed_time()) +
        " seconds.\n");
  }
#endif

  settings::alpha_omega_timer.stop();
  Output::instance()->write(
      " Total Runtime: " +
      std::to_string(settings::alpha_omega_timer.elapsed_time()) +
      " seconds.\n");

  return 0;
}
