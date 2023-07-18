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
