// $Source$
//----------------------------------------------------------------------
// PoleMatrix
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file PoleMatrix.cpp
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates 
 * for a given date
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\PoleMatrix.hpp"

/**
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates 
 * for a given date
 * 
 * @param xp Pole coordinte(xp,yp)
 * @param yp Pole coordinte(xp,yp)
 * @return Matrix& Pole matrix
 */

Matrix& PoleMatrix(double xp, double yp) {
    Matrix& polemat = R_y(-xp) * R_x(-yp);

    return polemat;
}