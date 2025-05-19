// $Header$
//----------------------------------------------------------------------
// gibbs
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file gibbs.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función gibbs.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _GIBBS_
#define _GIBBS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include <cmath>
#include <tuple>

tuple<Matrix&, double, double, double,string> gibbs(Matrix r1, Matrix r2, Matrix r3);

#endif