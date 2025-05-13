// $Header$
//----------------------------------------------------------------------
// NutAngles
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file NutAngles.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función NutAngles.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _NUTANGLES_
#define _NUTANGLES_

#include "..\include\SAT_Const.hpp"
#include "..\include\matrix.hpp"

#include <tuple>
#include <cmath>

tuple<double,double> NutAngles (double Mjd_TT);

#endif