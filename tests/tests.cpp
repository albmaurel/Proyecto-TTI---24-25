#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\timediff.hpp"
#include "..\include\sign_.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\IERS.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Frac.hpp"
#include "..\include\Legendre.hpp"



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
int m_extract_vector_01() {
    int c = 4;
    
    Matrix A(c);
    A(1) = 0; A(2) = 2; A(3) = 8; A(4) = 0;
    
    Matrix Expected(3);
    Expected(1) = 0;
    Expected(2) = 2;
    Expected(3) = 8;

    Matrix R = extract_vector(A, 1, 3);  
    
    _assert(m_equals(Expected, R, 1e-10));  
    
    return 0;
}

int m_union_vector_01() {
    
    Matrix A(3);  
    A(1) = 1;
    A(2) = 2;
    A(3) = 3;
    
    Matrix B(2);  
    B(1) = 4;
    B(2) = 5;

    Matrix Expected(5);
    Expected(1) = 1;
    Expected(2) = 2;
    Expected(3) = 3;
    Expected(4) = 4;
    Expected(5) = 5;

    Matrix R = union_vector(A, B);  
    
    _assert(m_equals(Expected, R, 1e-10));  
    
    return 0;
}

int m_extract_row_01() {
    int f = 3, c = 3;
    
    Matrix A(f, c);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix Expected(c);
    Expected(1) = 4;
    Expected(2) = 5;
    Expected(3) = 6;

    Matrix R = extract_row(A, 2);  

    _assert(m_equals(Expected, R, 1e-10));  

    return 0;
}

int m_extract_column_01() {
    int f = 3, c = 3;

    Matrix A(f, c);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

    Matrix Expected(f);
    Expected(1) = 2;
    Expected(2) = 5;
    Expected(3) = 8;

    Matrix R = extract_column(A, 2);  

    _assert(m_equals(Expected, R, 1e-10));  

    return 0;
}

int m_assign_row_01() {
    int f = 4, c = 3; 
    
    Matrix A(f, c);  
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;
    A(4,1) = 10; A(4,2) = 11; A(4,3) = 12;


    Matrix NewRow(3);  
    NewRow(1) = 13;
    NewRow(2) = 14;
    NewRow(3) = 15;

    Matrix R = assign_row(A, NewRow, 2);

	Matrix B(f, c);  
	B(1,1) = 1; B(1,2) = 2; B(1,3) = 3;
	B(2,1) = 13; B(2,2) = 14; B(2,3) = 15;
	B(3,1) = 7; B(3,2) = 8; B(3,3) = 9;
	B(4,1) = 10; B(4,2) = 11; B(4,3) = 12;
	
    _assert(m_equals(B, R, 1e-10));  

    return 0;
}

int m_assign_column_01() {
    int f = 4, c = 3; 
    
    Matrix A(f, c);  
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;
    A(4,1) = 10; A(4,2) = 11; A(4,3) = 12;


    Matrix NewColumn(4);  
    NewColumn(1) = 13;
    NewColumn(2) = 14;
    NewColumn(3) = 15;
	NewColumn(4) = 16;

    Matrix R = assign_column(A, NewColumn, 1);

	Matrix B(f, c);  
    B(1,1) = 13; B(1,2) = 2; B(1,3) = 3;
    B(2,1) = 14; B(2,2) = 5; B(2,3) = 6;
    B(3,1) = 15; B(3,2) = 8; B(3,3) = 9;
    B(4,1) = 16; B(4,2) = 11; B(4,3) = 12;
	
    _assert(m_equals(B, R, 1e-10));  

    return 0;
}

int R_x_01() {
    Matrix B = R_x(2.0);
	
    Matrix A(3, 3);  
    A(1,1) = 1.0000    ; A(1,2) = 0; A(1,3) =  0;
    A(2,1) = 0; A(2,2) = -0.4161; A(2,3) =  0.9093;
    A(3,1) = 0; A(3,2) = -0.9093 ; A(3,3) = -0.4161;

	
    _assert(m_equals(A, B, 1e-10));  

    return 0;
}

int R_y_01() {
	
    Matrix B = R_y(2.0);
	
    Matrix A(3, 3);  
    A(1,1) = -0.4161    ; A(1,2) = 0; A(1,3) =  -0.9093;
    A(2,1) =  0; A(2,2) = 1.0000; A(2,3) = 0;
    A(3,1) = 0.9093 ; A(3,2) = 0; A(3,3) = -0.4161;

	
    _assert(m_equals(A, B, 1e-10));  

    return 0;
}

int R_z_01() {
	
    Matrix B = R_z(2.0);
	
    Matrix A(3, 3);  
    A(1,1) = 0.4161    ; A(1,2) = 0.9093; A(1,3) =  0;
    A(2,1) = -0.9093; A(2,2) = -0.4161; A(2,3) = 0;
    A(3,1) = 0; A(3,2) = 0; A(3,3) = 1.0000;

	
    _assert(m_equals(A, B, 1e-10));  

    return 0;
}

int AzElPa_01() {
	
    Matrix A(3);  
    A(1) = 100;
    A(2) = 200;
    A(3) = 300;
	
	double Az, El;
    Matrix dAds(3), dEds(3);
	
	double Az_expected = 0.4636;
    double El_expected = 0.1334;
    Matrix dAds_expected(3), dEds_expected(3);

    dAds_expected(1) = 0.0004;
    dAds_expected(2) = -0.0002;
    dAds_expected(3) = 0.0;

    dEds_expected(1) = -0.0000264;
    dEds_expected(2) = -0.0000527;
    dEds_expected(3) = 0.0004393;
	
	AzElPa(A, Az, El, dAds, dEds);

	_assert(fabs(Az_expected - Az) < 1e-10);
    _assert(fabs(El_expected - El) < 1e-10);
	_assert(m_equals(dAds_expected, dAds, 1e-10));
	_assert(m_equals(dEds_expected, dEds, 1e-10));
	
    return 0;
}

int Mjday_01() {
    double mjd = Mjday(2024, 4, 17, 0, 0, 0); 
	double esperado = 60417.0;

	_assert(fabs(mjd - esperado) < 1e-10);
    return 0;
}

int Mjd_TDB_01() {
    double mjdtdb = Mjday_TDB(60417);
	double esperado = 60417.0;

	_assert(fabs(mjdtdb - esperado) < 1e-10);
    return 0;
}

int Position_01() {
    Matrix Pos = Position(0,1,1);
	
	Matrix A(3);  
    A(1) = 3454318.94409859;
    A(2) = 0;
    A(3) = 5343768.70088094;
	
	_assert(m_equals(A, Pos, 1e-10));

    return 0;
}

int timediff_01() {
    double UT1_UTC = -0.3341;
    double TAI_UTC = 37;

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    _assert(fabs(UT1_TAI + 37.3341) < 1e-10);
    _assert(fabs(UTC_GPS + 18.0) < 1e-10);
    _assert(fabs(UT1_GPS + 18.3341) < 1e-10);
    _assert(fabs(TT_UTC - 69.184) < 1e-10);
    _assert(fabs(GPS_UTC - 18.0) < 1e-10);
	
	return 0;

}

int sign__01(){
	double res=sign_(5,-3);
	double a =-5;
	_assert(fabs(res + a) < 1e-10);
	
	return 0;

}

int TimeUpdate_01(){
	Matrix P(2, 2);
    P(1, 1) = 1;
    P(1, 2) = 2;
    P(2, 1) = 3;
    P(2, 2) = 4;

    Matrix Phi(2, 2);
    Phi(1, 1) = 0.5;
    Phi(1, 2) = 0;
    Phi(2, 1) = 0;
    Phi(2, 2) = 0.5;

    double Qdt = 0.1;

    Matrix P_updated = TimeUpdate(P, Phi, Qdt);
	
	Matrix Res(2, 2);
    Phi(1, 1) = 0.25;
    Phi(1, 2) = 0.5;
    Phi(2, 1) = 0.75;
    Phi(2, 2) = 1.0;
	
	_assert(m_equals(Res, P_updated, 1e-10));
	
	return 0;
}

int AccelPointMass_01() {
    Matrix r(3);
    r(1) = 7000e3;
    r(2) = 0.0;
    r(3) = 0.0;

    Matrix s(3);
    s(1) = 384400e3;
    s(2) = 0.0;
    s(3) = 0.0;

    double GM = 4902.800066e9;

    Matrix a = AccelPointMass(r, s, GM);

    Matrix expected(3);
    expected(1) = 1.0e-05 *0.1242;
    expected(2) = 0.0;
    expected(3) = 0.0;

    _assert(m_equals(a, expected, 1e-12));

    return 0;
}

int Cheb3D_01() {
    double t = 0.5;
    int N = 3;
    double Ta = 0;
    double Tb = 1;

    Matrix Cx(3); Cx(1) = 1; Cx(2) = 2; Cx(3) = 3;
    Matrix Cy(3); Cy(1) = 0; Cy(2) = 1; Cy(3) = 0;
    Matrix Cz(3); Cz(1) = 5; Cz(2) = 0; Cz(3) = -1;

    Matrix ChebApp = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    Matrix expected(3);
    expected(1) = -2.0;   
    expected(2) = 0;
    expected(3) = 6;

    _assert(m_equals(ChebApp, expected, 1e-10));

    return 0;
}

int Frac_01() {
    _assert(fabs(Frac(5.75) - 0.7500) < 1e-10);
    return 0;
}

int MeanObliquity_01() {
    double Mjd_TT = 58000.0;
    double expected = 0.40905268985035; 

    double result = MeanObliquity(Mjd_TT);

    _assert(fabs(result - expected) < 1e-6);

    return 0;
}

int Legendre_01() {
    int n = 3;
    int m = 2;
    double fi = pi / 6.0;  

    Matrix pnm, dpnm;

    Legendre(n, m, fi, pnm, dpnm);

    Matrix expected_pnm(4, 4);  
    expected_pnm(1, 1) = 1.0;
	expected_pnm(2, 1) = 0.866025403784438;
    expected_pnm(2, 2) = 1.5;
    expected_pnm(3, 1) = -0.279508497187474;
    expected_pnm(3, 2) = 1.677050983124842;
	expected_pnm(3, 3) = 1.452368754827781;
	expected_pnm(4, 1) = -1.157516198590758;
    expected_pnm(4, 2) = 0.350780380010056;
    expected_pnm(4, 3) = 1.921303268617425;
	expected_pnm(4, 4) = 1.358566569955260;

    Matrix expected_dpnm(4, 4);  
    expected_pnm(1, 1) = 0;
	expected_pnm(2, 1) = 1.500000000000000;
    expected_pnm(2, 2) = -0.866025403784438;
    expected_pnm(3, 1) = 2.904737509655563;
    expected_pnm(3, 2) = 1.936491673103709;
	expected_pnm(3, 3) = -1.677050983124842;
	expected_pnm(4, 1) = 0.859232942804219;
    expected_pnm(4, 2) = 5.873171257932123;
    expected_pnm(4, 3) = 1.109264959331179;
	expected_pnm(4, 4) = -2.353106324627087;

	_assert(m_equals(pnm, expected_pnm, 1e-10));
	_assert(m_equals(dpnm, expected_dpnm, 1e-10));
		
	return 0;
}

int EccAnom_01() {
    double M = 1.0;  
    double e = 0.5;  

    double expected_E = 1.49870113351785;


    double calculated_E=EccAnom(M, e); // Adjusted to match the expected number of arguments

    _assert(fabs(expected_E - calculated_E) < 1e-6);

    return 0;
}

int NutAngles_01() {
    double Mj_dd = 55;  
    double dpsi = -1.168829531617437e-07;  

    double deps = -2.478350619864799e-08;


    double dpsi1, deps1;
    NutAngles(55, dpsi1, deps1);

    _assert(fabs(dpsi - deps1) < 1e-6);
    _assert(fabs(deps - dpsi1) < 1e-6);


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
	_verify(m_extract_vector_01);
	_verify(m_union_vector_01);  
	_verify(m_extract_row_01); 
    _verify(m_extract_column_01);
	_verify(m_assign_row_01);
	_verify(m_assign_column_01);
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(AzElPa_01);
	_verify(Mjday_01);
	_verify(Mjd_TDB_01);
	_verify(Position_01);
	_verify(timediff_01);
	_verify(sign__01);
	_verify(TimeUpdate_01);
	_verify(AccelPointMass_01);
	_verify(Cheb3D_01);
	_verify(Frac_01);
	_verify(MeanObliquity_01);
	_verify(EccAnom_01);
	_verify(NutAngles_01);


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
