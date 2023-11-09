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
#ifndef KAHAN_H
#define KAHAN_H


#include <simulation/particle.hpp>
/*
    This is actually the Kahan Summation Algorithm
*/
template<class InputIt, class T>
T kahan(InputIt first, InputIt last, T init) {
    T sum = init;
    T c = 0.;

    for (auto input = first; input != last; input++) {
        T y = *first - c;
        T t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

template<class InputIt, class T>
T kahanBabushka(InputIt first, InputIt last, T init) {
    T sum = init;
    T cs = 0.;
    T ccs = 0.;
    T c = 0.;
    T cc = 0.;

    for (auto input = first; input != last; input++) {
        T t = sum + *input;
        if(std::abs(sum) >= std::abs(*input))
            c = (sum - t) + *input;
        else
            c = (*input - t) + sum;
        sum = t;
        t = cs + c;
        if(std::abs(cs) >= std::abs(c))
            cc = (cs - t) + c;
        else
            cc = (c - t) + cs;
        cs = t;
        ccs = ccs + cc;
    }
    return sum + cs + ccs;
}

std::tuple<double,double,int,int> kahanBabushkaTuple(std::vector<BankedParticle>::iterator first, std::vector<BankedParticle>::iterator last) {

    double sum1 = 0.;
    double sum2 = 0.;

    int sum3 = 0;
    int sum4 = 0;

    double cs1 = 0.;
    double ccs1 = 0.;
    double c1 = 0.;
    double cc1 = 0.;

    double cs2 = 0.;
    double ccs2 = 0.;
    double c2 = 0.;
    double cc2 = 0.;

    for (auto input = first; input != last; input++) {
        if(input->wgt > 0.)
        {
            sum3++;
            double t = sum1 + input->wgt;

            if(std::abs(sum1) >= std::abs(input->wgt))
                c1 = (sum1 - t) + input->wgt;
            else
                c1 = (input->wgt - t) + sum1;

            sum1 = t;
            t = cs1 + c1;

            if(std::abs(cs1) >= std::abs(c1))
                cc1 = (cs1 - t) + c1;
            else
                cc1 = (c1 - t) + cs1;
            cs1 = t;
            ccs1 = ccs1 + cc1;
        }
        else
        {
            sum4++;
            double t = sum2 + input->wgt;

            if(std::abs(sum2) >= std::abs(input->wgt))
                c2 = (sum2 - t) + input->wgt;
            else
                c2 = (input->wgt - t) + sum2;
            sum2 = t;
            t = cs2 + c2;
            if(std::abs(cs2) >= std::abs(c2))
                cc2 = (cs2 - t) + c2;
            else
                cc2 = (c2 - t) + cs2;
            cs2 = t;
            ccs2 = ccs2 + cc2;
        }
    }

    double r1 = sum1 + cs1 + ccs1;
    double r2 = (sum2 + cs2 + ccs2) * -1.0;
    return std::tuple<double,double,int,int>(r1,r2,sum3,sum4);

}
#endif