// $Source$
//----------------------------------------------------------------------
// Frac
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Frac.cpp
 * @brief Computes the fractional part of a given number.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------

#include "..\include\Frac.hpp"


/**
 * @brief Computes the fractional part of a given number.
 * 
 * This function calculates the fractional part of a number by subtracting 
 * the largest integer less than or equal to the number from the number itself.
 * 
 * @param x The input number.
 * @return double The fractional part of the input number (x - floor(x)).
 */

double Frac(double x) {
    return x - floor(x);
}