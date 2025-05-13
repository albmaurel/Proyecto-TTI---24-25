// $Header$
//----------------------------------------------------------------------
// Position
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Position.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función Position.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef _POSITION_
#define _POSITION_
#include "..\include\SAT_Const.hpp"
#include "..\include\matrix.hpp"

#include <cmath>

Matrix& Position (double lon, double lat, double h);

#endif