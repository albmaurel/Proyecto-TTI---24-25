// $Header$
//----------------------------------------------------------------------
// Cheb3D
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Cheb3D.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función Cheb3D.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _CHEB3D_
#define _CHEB3D_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix& Cheb3D(double t, int N, double Ta, double Tb,  Matrix& Cx,  Matrix& Cy,  Matrix& Cz);

#endif