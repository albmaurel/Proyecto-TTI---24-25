// $Header$
//----------------------------------------------------------------------
// timediff
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file timediff.hpp
*	@brief Este archivo de cabecera contiene una implementación de la acción timediff.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef _TIMEDIFF_
#define _TIMEDIFF_
#include <cmath>

void timediff(double UT1_UTC, double TAI_UTC,
              double& UT1_TAI, double& UTC_GPS,
              double& UT1_GPS, double& TT_UTC,
              double& GPS_UTC);

#endif