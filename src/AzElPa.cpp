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
 * @brief Calculates the azimuth and elevation angles, along with their partial derivatives, from a position vector.
 *
 * This function computes the azimuth (Az) and elevation (El) angles based on the provided position vector `s`.
 * It also calculates the partial derivatives of Az and El with respect to the components of `s`.
 *
 * @param s Reference to a 3-element position vector (Matrix&).
 * @return A tuple containing:
 *   - Azimuth angle (double)
 *   - Elevation angle (double)
 *   - Reference to the partial derivatives of Azimuth with respect to `s` (Matrix&)
 *   - Reference to the partial derivatives of Elevation with respect to `s` (Matrix&)
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