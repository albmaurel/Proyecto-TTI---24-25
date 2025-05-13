// $Header$
//----------------------------------------------------------------------
// TimeUpdate
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file TimeUpdate.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función TimeUpdate.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef _TIME_UPDATE_
#define _TIME_UPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix TimeUpdate( Matrix& P,  Matrix& Phi, double Qdt = 0.0);

#endif