// $Source$
//----------------------------------------------------------------------
// sign_
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file sign_.cpp
 * @brief Returns the absolute value of a with the sign of b.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\sign_.hpp"


/**
 * @brief Returns the absolute value of a with the sign of b.
 * 
 * @param a The value whose absolute value is to be computed.
 * @param b The value whose sign is to be applied to the result.
 * @return double The absolute value of `a` with the sign of `b`.
 */

double sign_(double a, double b) {
    return (b > 0.0) ? fabs(a) : -fabs(a);
}