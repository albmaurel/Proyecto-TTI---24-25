// $Header$
//----------------------------------------------------------------------
// PrecMatrix
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file PrecMatrix.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función PrecMatrix.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\R_z.hpp"
#include "..\include\R_y.hpp"
#include "..\include\SAT_Const.hpp"

Matrix& PrecMatrix(double Mjd_1, double Mjd_2);

#endif