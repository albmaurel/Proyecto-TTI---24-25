// $Source$
//----------------------------------------------------------------------
// MeanObliquity
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file MeanObliquity.cpp
 * @brief Computes the mean obliquity of the ecliptic
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\MeanObliquity.hpp"


/**
 * @brief Computes the mean obliquity of the ecliptic
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double - Mean obliquity of the ecliptic [rad]
 */

double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - SAT_Const::MJD_J2000) / 36525.0;

    double arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    return SAT_Const::Rad * (arcsec / 3600.0);  
}