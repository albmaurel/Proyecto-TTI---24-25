// $Source$
//----------------------------------------------------------------------
// AccelPointMass
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file AccelPointMass.cpp
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\AccelPointMass.hpp"


/**
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @param r Satellite position vector
 * @param s Point mass position vector  
 * @param GM Gravitational coefficient of point mass
 * @return Matrix& Acceleration (a=d^2r/dt^2)
 */

Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM) {
    Matrix d = r - s;

    Matrix& a = d / (pow(norm(d), 3.0)) + s / (pow(norm(s), 3.0));
    a = a * (-GM);

    return a;
}