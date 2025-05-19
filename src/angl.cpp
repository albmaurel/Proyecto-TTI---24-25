// $Source$
//----------------------------------------------------------------------
// angl
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file angl.cpp
 * @brief This file contains the implementation of the angl function.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\angl.hpp"


/**
 * @brief Calculates the angle (in radians) between two vectors.
 * 
 * This function computes the angle between two vectors using the dot product formula.
 * If either vector has a near-zero magnitude, the function returns a predefined undefined value.
 * 
 * @param vec1 The first vector as a Matrix.
 * @param vec2 The second vector as a Matrix.
 * @return The angle in radians between vec1 and vec2, or a large undefined value if the angle cannot be determined.
 */
double angl(Matrix vec1, Matrix vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;
    double magv1 = norm(vec1);
    double magv2 = norm(vec2);
    double theta  = 0.0;
	double temp;
    if (magv1*magv2 > small^2) {
        temp= dot(vec1,vec2) / (magv1*magv2);
        if (abs( temp ) > 1.0) {
            double sign = temp > 0 ? 1 : -1;
            temp= sign* 1.0;
        }
         theta= acos( temp );
    } else {
        theta= undefined;
    }

    return theta;
}
