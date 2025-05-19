// $Source$
//----------------------------------------------------------------------
// unit
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file unit.cpp
 * @brief This file contains the implementation of the unit function.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\unit.hpp"

/**
 * @brief Calculates a unit vector given the original vector
 * 
 * If a zero vector is input, the vector is set to zero.
 * 
 * @param vec Vector to be converted to a unit vector
 * @return Matrix& Unit vector of the input vector 
 */

Matrix& unit(Matrix vec) {
    double small = 0.000001;
    double magv = norm(transpose(vec));
    Matrix* outvec = new Matrix(3);
    if ( magv > small ) {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= vec(i)/magv;
        }
    } else {
        for (int i=1; i <= 3; i++) {
            (*outvec)(i)= 0.0;
        }
    }

    return *outvec;

}