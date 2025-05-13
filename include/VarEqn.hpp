// $Header$
//----------------------------------------------------------------------
// VarEqn
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file VarEqn.hpp
*	@brief Este archivo de cabecera contiene una implementación de la función VarEqn.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//---------------------------------------------------------------------
#ifndef _VAREQN_
#define _VAREQN_

#include "../include/matrix.hpp"
#include "../include/global.hpp"
#include "../include/SAT_Const.hpp"
#include "../include/IERS.hpp"
#include "../include/timediff.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/G_AccelHarmonic.hpp"


Matrix& VarEqn(double x, Matrix yPhi);

#endif