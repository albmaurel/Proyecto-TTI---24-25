// $Header$
//----------------------------------------------------------------------
// AzElPa
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file AzElPa.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función AzElPa.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef __AZELPA__
#define _AZELPA_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>


tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s);

#endif