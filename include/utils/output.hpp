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
#ifndef OUTPUT_H
#define OUTPUT_H

#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <utils/header.hpp>
#include <utils/mpi.hpp>
#include <vector>

//============================================================================
// Output Class (Singleton)
class Output {
 public:
  ~Output() = default;

  static void set_output_filename(std::string fname);

  static std::shared_ptr<Output> instance();

  template <class T>
  void write(T message) {
    if (mpi::rank == 0) {
      write_mutex.lock();
      std::cout << message << std::flush;
      if (out.is_open()) out << message << std::flush;
      write_mutex.unlock();
    }
  }

  void save_warning(std::string message);

  void write_saved_warnings();

  void write_error(std::string message);

 private:
  Output();
  static std::shared_ptr<Output> output_instance;
  static std::string output_filename;
  static std::mutex instance_mutex;
  static std::mutex write_mutex;
  static std::mutex save_mutex;
  std::vector<std::string> warnings_for_latter;
  static bool fname_set;
  std::ofstream out;

};  // Output

//============================================================================
// Non-Member Functions for Output
void print_header();

std::string current_date_time();

#endif
