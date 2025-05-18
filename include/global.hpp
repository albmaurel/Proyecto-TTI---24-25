// $Header$
//----------------------------------------------------------------------
// global
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file global.hpp
*	@brief Este archivo de cabecera contiene una implementación de las acciones eop19620101,GGM03S,DE430Coeff,auxparam.
*	
*	@author Alberto Maurel Mendizábal
*	@bug No known bugs.
*/ 
//----------------------------------------------------------------------
#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include <string.h>
#include <cmath>

extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;
extern Param AuxParam;
extern Matrix Obs;

void eop19620101(int c);
void GGM03S();
void DE430Coeff();
void auxparam();
void GEOS3(int f);



#endif 