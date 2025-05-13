// $Header$
//----------------------------------------------------------------------
// JPL_Eph_DE430
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file JPL_Eph_DE430.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función JPL_Eph_DE430.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include "../include/matrix.hpp"
#include "../include/global.hpp"
#include "../include/Cheb3D.hpp"

#include <tuple>
tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);


#endif