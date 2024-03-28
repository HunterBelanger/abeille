#ifndef SYMMETRIC_MATRIX_INDEX_H
#define SYMMETRIC_MATRIX_INDEX_H

#include <iostream>

size_t symm_matrix_length(const size_t& order){
    return order * (order +1) * 0.5 ;
}

size_t symm_linear_index(const size_t& order, const size_t& i, const size_t& j){
    // i is for rows and j is for columns, and belongs to [0, order-1]
        if ( j < i ){
            return symm_linear_index(order, j, i);
        }else{
            return order * i - 0.5 * (i*i + i) + j ;
        }
}


#endif