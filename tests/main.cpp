#include "..\include\matrix.hpp"
#include <iostream>

using namespace std;

int main() {
	eop19620101(); // c=21413
	cout<<eopdata<<endl;
    Matrix M1(3, 2);
	M1(1,1) = 5;
	M1(1,2) = 3;
	
    Matrix M2(3, 2);
	M2(1,1) = -3;
	
    int f1 = 3, c1 = 3; 
    int f2 = 3, c2 = 3;  

    Matrix A(f1, c1);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 0; A(2,2) = 1; A(2,3) = 4;
    A(3,1) = 5; A(3,2) = 6; A(3,3) = 0;
	
    Matrix M3 = M1 - M2;

    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	
	cout << M1(1,1) << "\n";
	
	Matrix M4 = transpose(M1); 
    cout << "M1 Transposed:\n" << M4 << "\n";
	
	Matrix M5 = eye(4);
    cout << "Identity Matrix (3x3):\n" << M5 << "\n";

	try {
        Matrix M6 = inv(A); 
        cout << "Inverse of :A\n" << M6 << "\n";
    } catch (const std::exception &e) {
        cout << "Error: " << e.what() << "\n";  
    }
    return 0;
}