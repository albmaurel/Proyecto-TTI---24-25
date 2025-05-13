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

void AzElPa( Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {
    double rho = sqrt(s(1) * s(1) + s(2) * s(2));

    // Angles
    Az = atan2(s(1), s(2)); // Azimuth

    if (Az < 0.0) {
        Az += SAT_Const::pi2;
    }

    El = atan(s(3) / rho); // Elevation

    // Partials
    dAds(1) = s(2) / (rho * rho);
    dAds(2) = -s(1) / (rho * rho);
    dAds(3) = 0.0;

    dEds(1) = -s(1) * s(3) / rho;
    dEds(2) = -s(2) * s(3) / rho;
    dEds(3) = rho;

    // Normalize dEds by the dot product of s
    double dot_product = dot(s,s);
    dEds = dEds / dot_product;
}