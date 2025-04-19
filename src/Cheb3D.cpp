#include "..\include\matrix.hpp"
#include "..\include\Cheb3D.hpp"
#include <cmath>

Matrix& Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz){
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
	Matrix f3(c);
	f3(1) = Cx(i); f3(2) = Cy(i); f3(3) = Cz(i);
    f1 = f1*(2*tau)-f2+f3;
    f2 = old_f1;
	}
	
	Matrix f4(c);
	f4(1) = Cx(1); f4(2) = Cy(1); f4(3) = Cz(1);
	ChebApp = f1*tau-f2+f4;
	
	return ChebApp 
	
}