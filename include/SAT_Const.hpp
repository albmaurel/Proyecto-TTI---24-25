// $Header$
//----------------------------------------------------------------------
// SAT_Const
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file SAT_Const.hpp
*	@brief Este archivo de cabecera contiene una implementación del namespace SAT_Const.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef INITIAL_ORBIT_DETERMINATION_C___SAT_CONST_H
#define INITIAL_ORBIT_DETERMINATION_C___SAT_CONST_H

#ifndef M_PI
#define M_PI 3.14159265358979324
#endif
#include <cmath>

namespace SAT_Const {
    // Mathematical constants
    constexpr double pi2       = M_PI * 2.0;                 // 2π
  
    constexpr double Rad       = M_PI / 180.0;               // Radians per degree

    constexpr double Deg       = 180.0 / M_PI;               // Degrees per radian

    constexpr double Arcs      = 3600.0 * 180.0 / M_PI;        // Arcseconds per radian

    constexpr double MJD_J2000 = 51544.5;                  // Modified Julian Date of J2000

    constexpr double T_B1950   = -0.500002108;             // Epoch B1950

    constexpr double c_light   = 299792458.0;              // Speed of light [m/s]

    constexpr double AU        = 149597870700.0;           // Astronomical unit [m]

    constexpr double R_Earth   = 6378.1363e3;              // Earth's radius [m]

    constexpr double f_Earth   = 1.0 / 298.257223563;      // Flattening
 
    constexpr double R_Sun     = 696000e3;                 // Sun's radius [m]
 
    constexpr double R_Moon    = 1738e3;                   // Moon's radius [m]

    constexpr double omega_Earth = 15.04106717866910 / 3600 * Rad; // [rad/s]

    constexpr double GM_Earth    = 398600.435436e9;                           // [m^3/s^2]

    constexpr double GM_Sun      = 132712440041.939400e9;                     // [m^3/s^2]

    constexpr double GM_Moon     = GM_Earth / 81.30056907419062;             // [m^3/s^2]

    constexpr double GM_Mercury  = 22031.780000e9;                            // [m^3/s^2]

    constexpr double GM_Venus    = 324858.592000e9;                           // [m^3/s^2]

    constexpr double GM_Mars     = 42828.375214e9;                            // [m^3/s^2]

    constexpr double GM_Jupiter  = 126712764.800000e9;                        // [m^3/s^2]

    constexpr double GM_Saturn   = 37940585.200000e9;                         // [m^3/s^2]

    constexpr double GM_Uranus   = 5794548.600000e9;                          // [m^3/s^2]

    constexpr double GM_Neptune  = 6836527.100580e9;                          // [m^3/s^2]

    constexpr double GM_Pluto    = 977.0000000000009e9;                       // [m^3/s^2]

    constexpr double P_Sol       = 1367 / c_light;           // [N/m^2]
	
	constexpr double eps        = 2.22044604925031e-16;
}

#endif