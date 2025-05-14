// $Source$
//----------------------------------------------------------------------
// Position
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Position.cpp
 * @brief Calculate geocentric position from geodetic coordinates
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\Position.hpp"


/**
 * @brief Calculate geocentric position from geodetic coordinates
 * 
 * @param lon Longitude in radians
 * @param lat Latitude in radians  
 * @param h Altitude in meters
 * @return Matrix& Position vector in Earth-fixed frame
 */

Matrix& Position (double lon, double lat, double h){
	
	
	double R_equ = SAT_Const::R_Earth;
	double f     = SAT_Const::f_Earth;

	double e2     = f*(2.0-f);   // Square of eccentricity
	double CosLat = cos(lat);    // (Co)sine of geodetic latitude
	double SinLat = sin(lat);

	// Position vector 
	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
	Matrix& r = zeros(3);
	r(1) =  (N+h)*CosLat*cos(lon);
	r(2) =  (N+h)*CosLat*sin(lon);
	r(3) =  ((1.0-e2)*N+h)*SinLat;
	
	return r;
}