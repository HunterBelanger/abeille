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
#include <plotting/plotter.hpp>
#include <utils/constants.hpp>
#include <utils/header.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

//remove later
#include <simulation/approximate_mesh_cancelator.hpp>
#include <utils/position.hpp>
#include <unordered_map>
#include <ndarray.hpp>

#include <docopt.h>

#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#ifdef ABEILLE_USE_OMP
#include <omp.h>
#endif

void signal_handler(int signal) {
  if (!simulation->signaled) {
    // Signal has been recieved. Set the signaled
    // flag in simulation. It will stop after the
    // current generation's particles have all been
    // transported as we can't break inside an
    // OpenMP loop.
#ifdef ABEILLE_USE_OMP
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
#ifdef ABEILLE_USE_OMP
#pragma omp critical
#endif
    {
      Output::instance().write("\n !! >>> FORCE KILLED BY USER <<< !!\n");
      mpi::abort_mpi();
      std::exit(1);
    }
  }
}

bool exists(std::string fname) {
  std::ifstream file(fname);
  return file.good();
}


/*
TODO:

1. mpi.hpp add the short type so that vector version can use it
2. turn the vectors into NDArray
3. Revert things to how they were before testing 
4. Add the option for loop or vector


*/
int main(int argc, char** argv) {
  
  settings::alpha_omega_timer.start();

  mpi::initialize_mpi(&argc, &argv);

  std::atexit(mpi::finalize_mpi);

  Position start(-1.,-1.,-1.);
  Position end(1.,1.,1.);
  ApproximateMeshCancelator cancelator(start,end,2,1,1);
  std::vector<BankedParticle> v;
  std::cout << " made cancelator" << "\n";
    if(mpi::rank == 0)
    {
      v.emplace_back();
      v.back().r = Position(-0.5,0,0);
      v.back().wgt = 0.5;
      v.back().wgt2 = 0.5;
      v.emplace_back();
      v.back().r = Position(-0.3,0,0);
      v.back().wgt = 0.6;
      v.back().wgt2 = 0.5;
      v.emplace_back();
      v.back().r = Position(-0.4,0,0);
      v.back().wgt = 0.7;
      v.back().wgt2 = 0.5;
      v.emplace_back();
      v.back().r = Position(0.4,0,0);
      v.back().wgt = 0.6;
      v.back().wgt2 = 0.5;
      v.emplace_back();
      v.back().r = Position(0.3,0,0);
      v.back().wgt = 0.7;
      v.back().wgt2 = 0.5;
      v.emplace_back();
      v.back().r = Position(0.5,0,0);
      v.back().wgt = 0.8;
      v.back().wgt2 = -0.2;

      for(auto &p : v)
        cancelator.add_particle(p);
    }
    else if(mpi::rank == 1)
    {
      v.emplace_back();
      v.back().r = Position(-0.5,0,0);
      v.back().wgt = -0.3;
      v.back().wgt2 = -0.9;
      v.emplace_back();
      v.back().r = Position(-0.3,0,0);
      v.back().wgt = -0.3;
      v.back().wgt2 = -0.2;
      v.emplace_back();
      v.back().r = Position(-0.4,0,0);
      v.back().wgt = -0.1;
      v.back().wgt2 = -0.2;
      v.emplace_back();
      v.back().r = Position(0.4,0,0);
      v.back().wgt = -0.4;
      v.back().wgt2 = -0.5;
      v.emplace_back();
      v.back().r = Position(0.3,0,0);
      v.back().wgt = -0.9;
      v.back().wgt2 = -0.9;
      v.emplace_back();
      v.back().r = Position(0.5,0,0);
      v.back().wgt = -0.2;
      v.back().wgt2 = -0.1;
      
      for(auto &p : v)
        cancelator.add_particle(p);
    }
    /*
    else if(mpi::rank == 2)
    {
     
    }
    else if(mpi::rank == 3)
    {
     
    }
    */
  std::cout << " cancelating" << "\n";
  //cancelator.perform_cancellation_vector(settings::rng);
  cancelator.perform_cancellation_vector(settings::rng);
  std::cout << " cancelated" << "\n";

  if(mpi::rank == 0)
  {
    sleep(3);
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[0][0]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[0][0]->wgt2<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[1][0]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[1][0]->wgt2<< " \n";
    std::cout << " done " << "\n";
  }
if(mpi::rank == 1)
  {
    sleep(8);
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[0][0]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[0][0]->wgt2<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[1][0  ]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[1][0]->wgt2<< " \n";
    std::cout << " done " << "\n";
  }
  /*
if(mpi::rank == 2)
  {
    sleep(13);
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[2][0]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[2][0]->wgt2<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[2][1]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[2][1]->wgt2<< " \n";
    std::cout << " done " << "\n";

  }
if(mpi::rank == 3)
  {
    sleep(18);
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[3][0]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[3][0]->wgt2<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt: " << cancelator.bins[3][1]->wgt<< " \n";
    std::cout << " my rank is " << mpi::rank <<  " cancelator wgt2: " << cancelator.bins[3][1]->wgt2<< " \n";
    std::cout << " done " << "\n";
  }
  */
  
  /*
  // Make help message string for docopt
  std::string help_message = version_string + "\n" + info + "\n" + help;

  // Initialize docopt
  std::map<std::string, docopt::value> args =
      docopt::docopt(help_message, {argv + 1, argv + argc}, true,
                     version_string + "\n" + info);

  std::string output_filename;
  if (args["--output"]) {
    output_filename = args["--output"].asString();
    Output::set_output_filename(output_filename);
  }
  mpi::synchronize();

#ifdef ABEILLE_GUI_PLOT
  if (args["--gui-plot"].asBool() || args["--plot"].asBool()) {
#else
  if (args["--plot"].asBool()) {
#endif
    // If using OpenMP, get number of threads requested
#ifdef ABEILLE_USE_OMP
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
      // Set plotting mode to true. This way we don't open an HDF5 output file.
      settings::plotting_mode = true;

      print_header();
      Output::instance().write("\n");

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
        Output::instance().write(" FATAL ERROR: " + mssg + ".\n");
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
#ifdef ABEILLE_USE_OMP
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
    Output::instance().write("\n");

    // Parse input file
    bool parsed_file = false;
    try {
      parse_input_file(args["--input"].asString());
      parsed_file = true;
    } catch (const std::runtime_error& err) {
      std::string mssg = err.what();
      Output::instance().write(" FATAL ERROR: " + mssg + ".\n");
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
    Output::instance().write(
        " Total MPI time: " + std::to_string(mpi::timer.elapsed_time()) +
        " seconds.\n");

    if (mpi::rank == 0)
      Output::instance().h5().createAttribute("total-mpi-time",
                                              mpi::timer.elapsed_time());
  }
#endif

  settings::alpha_omega_timer.stop();
  Output::instance().write(
      " Total Runtime: " +
      std::to_string(settings::alpha_omega_timer.elapsed_time()) +
      " seconds.\n");

  if (mpi::rank == 0)
    Output::instance().h5().createAttribute(
        "total-run-time", settings::alpha_omega_timer.elapsed_time());

  return 0;
  */
}
