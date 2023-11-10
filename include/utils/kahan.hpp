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
    This is actually the Kahan-Babushka Summation Algorithm
*/
template<class InputIt, class T>
inline T kahan(InputIt first, InputIt last, T init) {
    T sum = init;
    T cs = 0.;
    T ccs = 0.;
    T c = 0.;
    T cc = 0.;

    for (auto input = first; input != last; input++) {
        auto input_i = *input;
        T t = sum + input_i;
        if(std::abs(sum) >= std::abs(input_i))
            c = (sum - t) + input_i;
        else
            c = (input_i - t) + sum;
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


// Returns tuple which contains two doubles and two ints
// doubles are positive and negative sum of banked particle weights in a vector 
// ints are the total counts of positive particles and negative particles
inline std::tuple<double,double,int,int> kahan_bank(std::vector<BankedParticle>::iterator first, std::vector<BankedParticle>::iterator last) {
    double sum_pos = 0.;
    double sum_neg = 0.;

    int count_pos = 0;
    int coung_neg = 0;

    double cs_pos = 0.;
    double ccs_pos = 0.;
    double c_pos = 0.;
    double cc_pos = 0.;

    double cs_neg = 0.;
    double ccs_neg = 0.;
    double c_neg = 0.;
    double cc_neg = 0.;

    for (auto input = first; input != last; input++) {
        auto input_i_wgt = input->wgt;
        if(input_i_wgt > 0.)
        {
            count_pos++;
            double t = sum_pos + input_i_wgt;

            if(std::abs(sum_pos) >= std::abs(input_i_wgt))
                c_pos = (sum_pos - t) + input_i_wgt;
            else
                c_pos = (input_i_wgt - t) + sum_pos;

            sum_pos = t;
            t = cs_pos + c_pos;

            if(std::abs(cs_pos) >= std::abs(c_pos))
                cc_pos = (cs_pos - t) + c_pos;
            else
                cc_pos = (c_pos - t) + cs_pos;
            cs_pos = t;
            ccs_pos = ccs_pos + cc_pos;
        }
        else
        {
            coung_neg++;
            double t = sum_neg + input_i_wgt;

            if(std::abs(sum_neg) >= std::abs(input_i_wgt))
                c_neg = (sum_neg - t) + input_i_wgt;
            else
                c_neg = (input_i_wgt - t) + sum_neg;
            sum_neg = t;
            t = cs_neg + c_neg;
            if(std::abs(cs_neg) >= std::abs(c_neg))
                cc_neg = (cs_neg - t) + c_neg;
            else
                cc_neg = (c_neg - t) + cs_neg;
            cs_neg = t;
            ccs_neg = ccs_neg + cc_neg;
        }
    }

    double r1 = sum_pos + cs_pos + ccs_pos;
    double r2 = (sum_neg + cs_neg + ccs_neg) * -1.0;
    return std::tuple<double,double,int,int>(r1,r2,count_pos,coung_neg);

}
#endif