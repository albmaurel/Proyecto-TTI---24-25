#include "..\include\matrix.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;
	
	//Starting value
    M = fmod(M, 2.0 * pi);

    double E;
    if (e < 0.8) {
        E = M;
    } else {
        E = pi;
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