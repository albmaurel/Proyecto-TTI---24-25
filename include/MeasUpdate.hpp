// $Header$
//----------------------------------------------------------------------
// MeasUpdate
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file MeasUpdate.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función MeasUpdate.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "..\include\matrix.hpp"
#include <tuple>

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n);

#endif