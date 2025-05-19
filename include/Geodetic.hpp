// $Header$
//----------------------------------------------------------------------
// Geodetic
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Geodetic.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función Geodetic.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _GEODETIC_
#define _GEODETIC_
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <tuple>

tuple<double, double, double> Geodetic(Matrix r);

#endif