// $Source$
//----------------------------------------------------------------------
// Legendre
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Legendre.cpp
 * @brief Computes the Legendre polynomials and their derivatives.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\Legendre.hpp"

/**
 * @brief Computes the Legendre polynomials and their derivatives.
 * 
 * This function calculates the associated Legendre polynomials (Pn,m) and their derivatives (dPn,m/dθ)
 * up to a given degree (n) and order (m) for a specified angle (fi) in radians.
 * 
 * @param n Maximum degree of the Legendre polynomials.
 * @param m Maximum order of the Legendre polynomials.
 * @param fi Angle in radians for which the polynomials are computed.
 * @return tuple<Matrix, Matrix> A tuple containing two matrices:
 *         - The first matrix contains the values of the Legendre polynomials.
 *         - The second matrix contains the values of their derivatives.
 */

tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi)
{
    Matrix& pnm = zeros(n+1, m+1);
    Matrix& dpnm = zeros(n+1, m+1);

    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2,2)=sqrt(3.0)*cosf(fi);
    dpnm(2,2)=-sqrt(3.0)*sinf(fi);
	int i=0;
    //diagonal coefficients
    for( i = 2; i<= n; i++){
        pnm(i+1,i+1)= sqrt((2.0*i+1.0)/(2*i))*cosf(fi)*pnm(i,i);
        dpnm(i+1,i+1)= sqrt((2.0*i+1.0)/(2*i))*((cosf(fi)*dpnm(i,i))-(sinf(fi)*pnm(i,i)));
    }
    // horizontal first step coefficients
    for( i = 1; i<= n; i++){ 
        pnm(i+1,i)= sqrt(2.0*i+1.0)*sinf(fi)*pnm(i,i);
    
        dpnm(i+1,i)= sqrt(2.0*i+1.0)*((cosf(fi)*pnm(i,i))+(sinf(fi)*dpnm(i,i)));
    }
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(1){
        for( i = k; i<= n; i++){         
            pnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sinf(fi)*pnm(i,j+1))
            -(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*pnm(i-1,j+1)));
            dpnm(i+1,j+1)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sinf(fi)*dpnm(i,j+1))
                +(sqrt(2*i-1.0)*cosf(fi)*pnm(i,j+1))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
    return tie(pnm, dpnm);
}