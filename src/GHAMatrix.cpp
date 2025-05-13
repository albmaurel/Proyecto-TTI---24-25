// $Source$
//----------------------------------------------------------------------
// GHAMatrix
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file GHAMatrix.cpp
 * @brief Transformation from true equator and equinox to Earth equator and 
 * Greenwich meridian system 
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\GHAMatrix.hpp"

/**
 * @brief Transformation from true equator and equinox to Earth equator and 
 * Greenwich meridian system 
 * 
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return Matrix& - Greenwich Hour Angle matrix
 */

Matrix& GHAMatrix(double Mjd_UT1) {
    
    Matrix& GHAmat = R_z( gast(Mjd_UT1) );

    return GHAmat;
}