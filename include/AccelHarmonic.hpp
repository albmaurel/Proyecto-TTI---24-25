// $Header$
//----------------------------------------------------------------------
// AccelHarmonic
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file AccelHarmonic.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función AccelHarmonic.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\global.hpp"
#include <cmath>

Matrix& AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif