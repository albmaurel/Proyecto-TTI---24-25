// $Header$
//----------------------------------------------------------------------
// G_AccelHarmonic
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file G_AccelHarmonic.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función G_AccelHarmonic.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _G_ACCELHARMONIC_
#define _G_ACCELHARMONIC_

#include "../include/matrix.hpp"
#include "../include/AccelHarmonic.hpp"

Matrix& G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max);


#endif