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
#ifndef RNG_H
#define RNG_H

#include <pcg_random.hpp>
#include <random>
#include <utils/constants.hpp>

class RNG {
 public:
  // -----------------------------------------------------------------------
  // rand()
  //   This advances the pcg32 generator and produces a double
  //   over the interval [0,1).
  // -----------------------------------------------------------------------
  static double rand(pcg32& rng) { return unit_dist(rng); }

  // -----------------------------------------------------------------------
  // uniform(double a, double b)
  //   Returns a random double with a uniform distribution over the
  //   interval [a, b)
  //
  //   f(x|a,b) = 1 / (b - a)
  // -----------------------------------------------------------------------
  static double uniform(pcg32& rng, double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // normal(double mu, double sigma)
  //   Returns a random double from the normal (Gaussian) distribution,
  //   defined by average value mu, and std sigma
  //
  //   f(x|mu,sigma) = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x - mu)/sigma)^2)
  // -----------------------------------------------------------------------
  static double normal(pcg32& rng, double mu, double sigma) {
    std::normal_distribution<double> dist(mu, sigma);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // exponential(double lambda)
  //   Returns a random double from the exponential distribution,
  //   defined by the average value 1/lambda.
  //
  //   f(x|lambda) = lambda * exp(-lambda * x)
  // -----------------------------------------------------------------------
  static double exponential(pcg32& rng, double lambda) {
    if (lambda == 0.) return INF;

    std::exponential_distribution<double> dist(lambda);
    return dist(rng);
  }

  // -----------------------------------------------------------------------
  // discrete(std::vector<double> weights)
  //   Returns an integer over range [0,weights.size() - 1]
  //   where probability of each integer is defined by the weight
  //
  //   P(i|w_0, w_1, ... w_k) = w_i / Sum[j = 1 to k](w_j)
  // -----------------------------------------------------------------------
  static int discrete(pcg32& rng, const double* begin, const double* end) {
    std::discrete_distribution<int> dist(begin, end);
    return dist(rng);
  }

  static int discrete(pcg32& rng, const std::vector<double>& weights) {
    std::discrete_distribution<int> dist(weights.begin(), weights.end());
    return dist(rng);
  }

 private:
  static std::uniform_real_distribution<double> unit_dist;
};  // RNG

#endif
