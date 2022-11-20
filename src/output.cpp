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
#include <string>
#include <utils/constants.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>

#ifdef _OPENMP
#include <omp.h>  // For print_header function
#endif

//============================================================================
// Initialization of static members of Output singleton
std::shared_ptr<Output> Output::output_instance = nullptr;
std::string Output::output_filename = "output.txt";
std::mutex Output::instance_mutex;
std::mutex Output::write_mutex;
std::mutex Output::save_mutex;
bool Output::fname_set = false;

//============================================================================
// Output Singleton Methods
Output::Output() : warnings_for_latter(), out() {
  if (mpi::rank == 0) out.open(output_filename);
}

void Output::set_output_filename(std::string fname) {
  instance_mutex.lock();
  if (fname_set == false) {
    output_filename = fname;
    fname_set = true;
  }
  instance_mutex.unlock();
}

std::shared_ptr<Output> Output::instance() {
  if (output_instance == nullptr) {
    instance_mutex.lock();
    if (output_instance == nullptr) {
      output_instance = std::shared_ptr<Output>(new Output());
    }
    instance_mutex.unlock();
  }
  return output_instance;
}

void Output::save_warning(std::string message) {
  save_mutex.lock();
  warnings_for_latter.push_back(message);
  save_mutex.unlock();
}

void Output::write_saved_warnings() {
  if (!warnings_for_latter.empty()) {
    std::cout << "\n";
    if (out.is_open()) out << "\n";
  }

  for (const auto& wrn : warnings_for_latter) {
    std::cout << " WARNING: " << wrn << "\n";
    if (out.is_open()) out << " WARNING: " << wrn << "\n";
  }

  if (!warnings_for_latter.empty()) {
    std::cout << "\n";
    if (out.is_open()) out << "\n";
  }

  warnings_for_latter.clear();
}

void Output::write_error(std::string message) {
  write_mutex.lock();
  std::cerr << message << std::flush;
  if (out.is_open()) out << message << std::flush;
  write_mutex.unlock();
}

//============================================================================
// Non-Member Functions
void print_header() {
  std::shared_ptr<Output> output = Output::instance();
  output->write(logo);
  output->write(header);

  // Version info
  std::string info = "";
  info += " Version             : " + std::string(ABEILLE_VERSION_STRING);
  if (DEVELOPMENT_VERSION) {
    info += " (Development)\n";
  } else
    info += "\n";

  // Build date
  info += " Build Date          : " + std::string(__DATE__) + " ";
  info += std::string(__TIME__) + "\n";

  // Current date and time
  info += " Date/Time           : " + current_date_time() + "\n";

#ifdef _OPENMP
  info +=
      " OpenMP Threads      : " + std::to_string(omp_get_max_threads()) + "\n";
#endif

#ifdef ABEILLE_USE_MPI
  info += " MPI Ranks           : " + std::to_string(mpi::size) + "\n";
#endif

  info += " Abeille Git Hash   : " ABEILLE_GIT_HASH "\n";

  output->write(info);
}

std::string current_date_time() {
  std::time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%b %e %Y %H:%M:%S", &tstruct);
  return buf;
}
