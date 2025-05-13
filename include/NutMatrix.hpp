// $Header$
//----------------------------------------------------------------------
// NutMatrix
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file NutMatrix.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función NutMatrix.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"

Matrix& NutMatrix(double Mjd_TT);

#endif