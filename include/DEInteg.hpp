// $Header$
//----------------------------------------------------------------------
// DEInteg
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file DEInteg.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función DEInteg.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _DEINTEG_
#define _DEINTEG_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"

#include <cmath>

Matrix& DEInteg(Matrix& f(double t, Matrix y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif