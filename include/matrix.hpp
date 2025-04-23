#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;
	Matrix();
    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
	Matrix(const int v_size);
	// Member operators
	double& operator () (const int v_size);
	double& operator () (const int row, const int column);
	Matrix& operator = (Matrix& matrix2);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix& matrix2);
	Matrix& operator / (Matrix& matrix2);
	Matrix& operator + (double scalar);
	Matrix& operator - (double scalar);
	Matrix& operator * (double scalar);
	Matrix& operator / (double scalar);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);
Matrix& zeros(const int n);
Matrix& eye(int size);
Matrix& transpose(Matrix &m);
Matrix& inv(Matrix &m);
double norm(Matrix& a);
double dot(Matrix& a,Matrix& b);
Matrix& cross(Matrix& a,Matrix& b);
Matrix& extract_vector(Matrix& m, const int n1, const int n2);
Matrix& union_vector(Matrix& m1, Matrix& m2);
Matrix& extract_column(Matrix& m, const int column);
Matrix& extract_row(Matrix& m, const int row);
Matrix& assign_column(Matrix& m, Matrix& v, const int column);
Matrix& assign_row(Matrix& m, Matrix& v, const int row);

#endif