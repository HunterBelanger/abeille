#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H

#include <iostream>
#include <vector>

// upper triangular matrix for standard deviation in each bin for FET
// upper triangular matrix is made with linear index to save the memory.
class UpperTriangularMatrix{
    public:
        UpperTriangularMatrix() : order(0), matrix(0, 0.0) {};

        UpperTriangularMatrix(int M): order(M), matrix(M*(M+1)/2, 0.0) {};

        ~UpperTriangularMatrix() = default;

        // Get the order of the upper triangular matrix
        int get_matrix_order()const { return order; }

        // set the particular matrix element to a value
        void set_element(const unsigned int i, const unsigned int j, double x) {  matrix[linear_index(i,j)] = x;  }

        // get the particular matrix element
        double get_element(const unsigned int i,const unsigned int j)const {  return matrix[linear_index(i,j)];   }

        // add a value to a certain element
        void add_into_elementconst(const unsigned int i, const unsigned int j, const double x) {  matrix[linear_index(i,j)] += x; }

        // get the linear index of the element
        unsigned int linear_index(unsigned int i, unsigned int j)const {
        // i is for rows and j is for columns, and belongs to [0, order-1]
            if ( j < i )
                return linear_index(j,i);
            else{
                if ( i < 0 || i >= order || j < 0 || j >= order)
                    std::cout<<"Fatal Error, The given i and j are incorrect.\n";    // put the error option to be sure
                //i--;
                //return order * (i+1) - ((i*i + i) / 2) + (j - (i+1));
                //return order * (i+1) - 0.5 * (i*i + i) - (order-j) ;
                return order * i - 0.5 * (i*i + i) + j ;
            }
        }

    public:
        int order;
        std::vector<double> matrix;
};


#endif