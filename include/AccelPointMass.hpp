// $Header$
//----------------------------------------------------------------------
// AccelPointMass
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file AccelPointMass.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función AccelPointMass.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM);

#endif