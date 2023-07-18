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
#ifndef PIXEL_H
#define PIXEL_H

#include <cstdint>
#include <iostream>

namespace plotter {

//======================================================================
// Pixel class to contain info for a pixel or color definition
class Pixel {
 public:
  Pixel(uint8_t i_r = 0xff, uint8_t i_g = 0xff, uint8_t i_b = 0xff)
      : r(i_r), g(i_g), b(i_b) {}
  ~Pixel() = default;

  uint8_t R() const { return r; }
  uint8_t G() const { return g; }
  uint8_t B() const { return b; }

  void set_R(uint8_t i_r) { r = i_r; }
  void set_G(uint8_t i_g) { g = i_g; }
  void set_B(uint8_t i_b) { b = i_b; }

 private:
  uint8_t r, g, b;
};  // Pixel

//========================================================================
// Non-Member Function Declarations
inline std::ostream& operator<<(std::ostream& out, const Pixel& p) {
  out << p.R() << p.G() << p.B();
  return out;
}

}  // namespace plotter

#endif
