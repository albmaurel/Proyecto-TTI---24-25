// $Source$
//----------------------------------------------------------------------
// AzElPa
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file AzElPa.cpp
 * @brief Computes azimuth, elevation, and their partial derivatives from local tangent coordinates.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\AzElPa.hpp"


/**
 * @brief Computes azimuth, elevation, and their partial derivatives from local tangent coordinates.
 *
 * This method calculates the azimuth (Az) and elevation (El) angles in radians, as well as their partial derivatives
 * with respect to the topocentric local tangent coordinates (East-North-Zenith frame).
 *
 * @param s Input matrix representing topocentric local tangent coordinates (East-North-Zenith frame).
 * @param Az Output double representing the azimuth angle in radians.
 * @param El Output double representing the elevation angle in radians.
 * @param dAds Output matrix representing the partial derivatives of azimuth with respect to s.
 * @param dEds Output matrix representing the partial derivatives of elevation with respect to s.
 */

tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s) {

    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    double Az = atan2(s(1),s(2));

    if (Az<0.0){ 
        Az = Az+SAT_Const::pi2;
    }

    double El = atan ( s(3) / rho );

    // Partials
    Matrix &dAds = zeros(3);
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;
	
    Matrix &dEds = zeros(3);
    dEds(1) = -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    dEds(3) = rho;
    dEds = dEds / dot(s,s);

    return tie(Az, El, dAds, dEds);
    
}