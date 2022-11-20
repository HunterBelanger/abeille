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
#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <ratio>

class Timer {
 public:
  using TimePoint = std::chrono::high_resolution_clock::time_point;

  Timer() {}

  void start() {
    if (!is_ticking_) {
      start_time = std::chrono::high_resolution_clock::now();
      is_ticking_ = true;
    }
  }

  void stop() {
    if (is_ticking_) {
      TimePoint stop_time = std::chrono::high_resolution_clock::now();
      is_ticking_ = false;

      double new_time = get_seconds(start_time, stop_time);
      previous_elapsed_time += new_time;
    }
  }

  void reset() {
    is_ticking_ = false;
    previous_elapsed_time = 0.;
  }

  double elapsed_time() const {
    TimePoint now = std::chrono::high_resolution_clock::now();
    double time_to_add = 0.;

    if (is_ticking_) {
      time_to_add = get_seconds(start_time, now);
    }

    return time_to_add + previous_elapsed_time;
  }

  bool is_ticking() const { return is_ticking_; }

 private:
  TimePoint start_time = std::chrono::high_resolution_clock::now();
  bool is_ticking_ = false;
  double previous_elapsed_time = 0.;

  double get_seconds(TimePoint t0, TimePoint t1) const {
    std::chrono::duration<double, std::nano> diff = t1 - t0;
    diff = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
    return diff.count() * 1.E-9;
  }
};

#endif
