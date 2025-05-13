// $Source$
//----------------------------------------------------------------------
// R_y
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file R_y.cpp
 * @brief Generates a 3x3 rotation matrix for a given angle around the y-axis.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\R_y.hpp"


/**
 * @brief Generates a 3x3 rotation matrix for a given angle around the y-axis.
 * 
 * 
 * @param angle The angle of rotation in radians.
 * @return Matrix& A reference to the resulting 3x3 rotation matrix.
 */

Matrix& R_y(double angle){
	double C = cos(angle);
	double S = sin(angle);
	Matrix& rotmat = zeros(3,3);

	rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
	rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
	rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;
	
	return rotmat;
}