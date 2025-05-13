// $Source$
//----------------------------------------------------------------------
// EccAnom
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file EccAnom.cpp
 * @brief Computes the eccentric anomaly for elliptic orbits
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\EccAnom.hpp"
/**
 * @brief Computes the eccentric anomaly for elliptic orbits
 *
 * @param M Mean anomaly in [rad]
 * @param e Eccentricity of the orbit [0,1]
 * @return Eccentric anomaly in [rad]
 */

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;
	
	//Starting value
    M = fmod(M, 2.0 * M_PI );

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI ;
    }

    double f = E - e * sin(E) - M;
    
    E = E - f / (1.0 - e * cos(E));
	
    double eps = 1e-10;
	//Iteration
    while (fabs(f) > 1e-2 * eps) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;

        if (i == maxit) {
			cout << "Convergence problems in EccAnom\n";
			exit(EXIT_FAILURE);
        }
    }

    return E;  
}