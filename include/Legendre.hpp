// $Header$
//----------------------------------------------------------------------
// Legendre
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Legendre.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función Legendre.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _LEGENDRE_
#define _LEGENDRE_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <iostream>

#include <cmath>
#include <tuple>



tuple<Matrix, Matrix> Legendre(int n, int m, double fi);

#endif