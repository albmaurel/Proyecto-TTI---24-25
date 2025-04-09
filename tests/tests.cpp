#include "..\include\matrix.h"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_zeros_02() {
    int c = 4;
	
	Matrix A(c);
	A(1) = 0; A(2) = 0; A(3) = 0; A(4) = 0;
	
	Matrix B = zeros(4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_mul_01() {
    int f1 = 3, c1 = 2;  
    int f2 = 2, c2 = 4; 

    Matrix A(f1, c1);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    A(3,1) = 5; A(3,2) = 6;

    Matrix B(f2, c2);
    B(1,1) = 1; B(1,2) = 0; B(1,3) = 2; B(1,4) = 3;
    B(2,1) = 4; B(2,2) = 5; B(2,3) = 6; B(2,4) = 7;

    Matrix C(f1, c2);
    C(1,1) = 9;  C(1,2) = 10;  C(1,3) = 14;  C(1,4) = 17;
    C(2,1) = 19; C(2,2) = 20;  C(2,3) = 30;  C(2,4) = 37;
    C(3,1) = 29; C(3,2) = 30;  C(3,3) = 46;  C(3,4) = 57;

    Matrix R = A * B;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_div_01() {
    int f1 = 3, c1 = 3; 
    int f2 = 3, c2 = 3;  

    Matrix A(f1, c1);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 0; A(2,2) = 1; A(2,3) = 4;
    A(3,1) = 5; A(3,2) = 6; A(3,3) = 0;

    Matrix B(f2, c2);
    B(1,1) = 1; B(1,2) = 0; B(1,3) = 0;
    B(2,1) = 0; B(2,2) = 1; B(2,3) = 0;
    B(3,1) = 0; B(3,2) = 0; B(3,3) = 1;

    Matrix C(f1, c2);
    C(1,1) = 1; C(1,2) = 2; C(1,3) = 3;
    C(2,1) = 0; C(2,2) = 1; C(2,3) = 4;
    C(3,1) = 5; C(3,2) = 6; C(3,3) = 0;

    Matrix R = A / B;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_eye_01() {
    int f = 3;
    Matrix A = eye(f); 

    Matrix C(f, f); 
    C(1,1) = 1; C(1,2) = 0; C(1,3) = 0;
    C(2,1) = 0; C(2,2) = 1; C(2,3) = 0;
    C(3,1) = 0; C(3,2) = 0; C(3,3) = 1;

    _assert(m_equals(A, C, 1e-10)); 

    return 0;
}

int m_transpose_01() {
    int f = 3, c = 2;
    Matrix A(f, c);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    A(3,1) = 5; A(3,2) = 6;

    Matrix C(c, f);
    C(1,1) = 1; C(1,2) = 3; C(1,3) = 5;
    C(2,1) = 2; C(2,2) = 4; C(2,3) = 6;

    Matrix B = transpose(A); 

    _assert(m_equals(B, C, 1e-10)); 

    return 0;
}

int m_inv_01() {
    int f = 3, c = 3;
    Matrix A(f, c);
    A(1, 1) = 1; A(1, 2) = 2; A(1, 3) = 3;
    A(2, 1) = 0; A(2, 2) = 1; A(2, 3) = 4;
    A(3, 1) = 5; A(3, 2) = 6; A(3, 3) = 0;

    Matrix A_inv(f, c);
    A_inv(1, 1) = -24; A_inv(1, 2) = 18; A_inv(1, 3) = 5;
    A_inv(2, 1) = 20; A_inv(2, 2) = -15; A_inv(2, 3) = -4;
    A_inv(3, 1) = -5; A_inv(3, 2) = 4; A_inv(3, 3) = 1;

    Matrix R = inv(A);

    _assert(m_equals(A_inv, R, 1e-10));

    return 0;
}

int m_add_scalar_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    double scalar = 5.0;

    Matrix C(f, c);
    C(1,1) = 6; C(1,2) = 7; C(1,3) = 8;
    C(2,1) = 9; C(2,2) = 10; C(2,3) = 11;
    C(3,1) = 12; C(3,2) = 13; C(3,3) = 14;

    Matrix R = A + scalar;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_sub_scalar_01() {
    int f = 3;
    int c = 3;

    Matrix A(f, c);
    A(1,1) = 10; A(1,2) = 20; A(1,3) = 30;
    A(2,1) = 40; A(2,2) = 50; A(2,3) = 60;
    A(3,1) = 70; A(3,2) = 80; A(3,3) = 90;

    double scalar = 5.0;

    Matrix C(f, c);
    C(1,1) = 5; C(1,2) = 15; C(1,3) = 25;
    C(2,1) = 35; C(2,2) = 45; C(2,3) = 55;
    C(3,1) = 65; C(3,2) = 75; C(3,3) = 85;

    Matrix R = A - scalar;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_mul_scalar_01() {
    int f = 2;
    int c = 2;

    Matrix A(f, c);
    A(1,1) = 2; A(1,2) = 3;
    A(2,1) = 4; A(2,2) = 5;

    double scalar = 3.0;

    Matrix C(f, c);
    C(1,1) = 6; C(1,2) = 9;
    C(2,1) = 12; C(2,2) = 15;

    Matrix R = A * scalar;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_div_scalar_01() {
    int f = 2;
    int c = 2;

    Matrix A(f, c);
    A(1,1) = 8; A(1,2) = 4;
    A(2,1) = 6; A(2,2) = 2;

    double scalar = 2.0;

    Matrix C(f, c);
    C(1,1) = 4; C(1,2) = 2;
    C(2,1) = 3; C(2,2) = 1;

    Matrix R = A / scalar;

    _assert(m_equals(C, R, 1e-10));

    return 0;
}

int m_norm_01() {
    Matrix a(3); 
    a(1) = 3; a(2) = 4; a(3) = 0;

    double expected = 5.0; 
    double result = norm(a);

    _assert(fabs(result - expected) < 1e-10);
    return 0;
}

int m_dot_01() {
    Matrix a(1, 3); 
    a(1,1) = 1; a(1,2) = 2; a(1,3) = 3;

    Matrix b(1, 3);
    b(1,1) = 4; b(1,2) = -5; b(1,3) = 6;

    double expected = 1*4 + 2*(-5) + 3*6; 
    double result = dot(a, b);

    _assert(fabs(result - expected) < 1e-10);
    return 0;
}

int m_cross_01() {
    Matrix a(1, 3);
    a(1,1) = 1; a(1,2) = 2; a(1,3) = 3;

    Matrix b(1, 3);
    b(1,1) = 4; b(1,2) = 5; b(1,3) = 6;

    Matrix expected(3);
    expected(1) = -3; 
    expected(2) = 6;  
    expected(3) = -3; 

    Matrix result = cross(a, b);

    _assert(m_equals(result, expected, 1e-10));
    return 0;
}

int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_zeros_01);
	_verify(m_zeros_02);
	_verify(m_mul_01);
	_verify(m_div_01);
	_verify(m_transpose_01);
	_verify(m_eye_01);
	_verify(m_inv_01);
	_verify(m_add_scalar_01);
    _verify(m_sub_scalar_01);
    _verify(m_mul_scalar_01);
    _verify(m_div_scalar_01);
	_verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
