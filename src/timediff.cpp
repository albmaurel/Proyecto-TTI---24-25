// $Source$
//----------------------------------------------------------------------
// timediff
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file timediff.cpp
 * @brief Computes various time differences between different time scales.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\timediff.hpp"

/**
 * @brief Computes various time differences between different time scales.
 * 
 * @param UT1_UTC Difference between UT1 and UTC time scales [s].
 * @param TAI_UTC Difference between TAI and UTC time scales [s].
 * @param UT1_TAI Output: Difference between UT1 and TAI time scales [s].
 * @param UTC_GPS Output: Difference between UTC and GPS time scales [s].
 * @param UT1_GPS Output: Difference between UT1 and GPS time scales [s].
 * @param TT_UTC Output: Difference between TT and UTC time scales [s].
 * @param GPS_UTC Output: Difference between GPS and UTC time scales [s].
 */
void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS,
              double& UT1_GPS, double& TT_UTC,
              double& GPS_UTC) {

    const double TT_TAI  = 32.184;          // TT-TAI time difference [s]
    const double GPS_TAI = -19.0;            // GPS-TAI time difference [s]
    const double TT_GPS  =  TT_TAI - GPS_TAI; // TT-GPS time difference [s]
    const double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

    UT1_TAI = UT1_UTC - TAI_UTC;  // UT1-TAI time difference [s]
    double UTC_TAI = -TAI_UTC;    // UTC-TAI time difference [s]

    UTC_GPS = UTC_TAI - GPS_TAI;  // UTC_GPS time difference [s]
    UT1_GPS = UT1_TAI - GPS_TAI;  // UT1-GPS time difference [s]

    TT_UTC  = TT_TAI - UTC_TAI;   // TT-UTC time difference [s]
    GPS_UTC = GPS_TAI - UTC_TAI;  // GPS-UTC time difference [s]
}