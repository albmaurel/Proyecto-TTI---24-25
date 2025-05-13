// $Header$
//----------------------------------------------------------------------
// IERS
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file IERS.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función IERS.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"

#include <cmath>
#include <tuple>
tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp);

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC);

#endif