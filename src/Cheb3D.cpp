// $Source$
//----------------------------------------------------------------------
// Cheb3D
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Cheb3D.cpp
 * @brief Chebyshev approximation of 3-dimensional vectors
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\Cheb3D.hpp"


/**
 * @brief Chebyshev approximation of 3-dimensional vectors
 *
 * @param t  Time value to evaluate at
 * @param N  Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval 
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate)
 * @return Matrix& Position vector at time t
 */

Matrix& Cheb3D(double t, int N, double Ta, double Tb,  Matrix& Cx,  Matrix& Cy,  Matrix& Cz){
    //Check validity
	if (t < Ta || t > Tb) {
		cout << "ERROR: Time out of range in Cheb3D::Value\n";
        exit(EXIT_FAILURE);
	}
	//Clenshaw algorithm
	double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);
	
	Matrix f1=zeros(3);
	Matrix f2=zeros(3);
	Matrix old_f1=zeros(3);
	
    for (int i = N; i >= 2; --i) {
    old_f1 = f1;
	Matrix f3(3);
	f3(1) = Cx(i); 
	f3(2) = Cy(i); 
	f3(3) = Cz(i);
	f1 = f1*(2*tau)-f2+f3;
	f2 = old_f1;
	}
	Matrix f4(3);
	f4(1) = Cx(1); 
	f4(2) = Cy(1); 
	f4(3) = Cz(1);
	Matrix& ChebApp = f1 * tau - f2 + f4;
	
	return ChebApp;
	
}