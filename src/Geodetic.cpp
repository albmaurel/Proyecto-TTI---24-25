// $Source$
//----------------------------------------------------------------------
// Geodetic
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Geodetic.cpp
 * @brief This file contains the implementation of the Geodetic function.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\Geodetic.hpp"


/**
 * @brief Converts Cartesian coordinates to geodetic coordinates (longitude, latitude, altitude).
 *
 * This function transforms a position vector in Cartesian coordinates (ECEF) to geodetic coordinates
 * (longitude, latitude, and altitude) based on the WGS-84 ellipsoid parameters.
 * The conversion uses an iterative method to account for the Earth's flattening.
 *
 * @param r A 3-element Matrix representing the Cartesian coordinates (X, Y, Z) in meters.
 * @return A tuple containing:
 *         - longitude (radians)
 *         - latitude (radians)
 *         - altitude above the ellipsoid (meters)
 */
tuple<double, double, double> Geodetic(Matrix r) {
    
    double R_equ = SAT_Const::R_Earth;
    double f     = SAT_Const::f_Earth;

    double epsRequ = 2.2204e-16*R_equ;        // Convergence criterion
    double e2      = f*(2.0-f);        // Square of eccentricity

    double X = r(1);                   // Cartesian coordinates
    double Y = r(2);
    double Z = r(3);
    double rho2 = X*X + Y*Y;           // Square of distance from z-axis
    double lon, lat, h;              // Longitude, latitude, altitude
    // Check validity of input data
    if (norm(r)==0.0) {
        printf ( " invalid input in Geodetic constructor\n" );
        
        lon = 0.0;
        lat = 0.0;
        h   = -R_equ;
    }

    // Iteration 
    double dZ = e2*Z;
    double ZdZ,Nh,N;

    while(1) {
         ZdZ    =  Z + dZ;
         Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        double SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
         N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        double dZ_new =  N*e2*SinPhi;
        if ( abs(dZ-dZ_new) < epsRequ ) {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;

    return tie(lon, lat, h);

}