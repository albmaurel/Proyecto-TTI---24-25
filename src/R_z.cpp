// $Source$
//----------------------------------------------------------------------
// R_z
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file R_z.cpp
 * @brief Generates a 3x3 rotation matrix for a given angle around the z-axis.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------

#include "..\include\R_z.hpp"

/**
 * @brief Generates a 3x3 rotation matrix for a given angle around the z-axis.
 * 
 * 
 * @param angle The angle of rotation in radians.
 * @return Matrix& A reference to the resulting 3x3 rotation matrix.
 */
Matrix& R_z(double angle){
	double C = cos(angle);
	double S = sin(angle);
	Matrix& rotmat = zeros(3,3);

	rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
	rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
	rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;
	
	return rotmat;
}