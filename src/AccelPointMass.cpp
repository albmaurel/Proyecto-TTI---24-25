#include "..\include\matrix.hpp"
#include "..\include\AccelPointMass.hpp"
#include <cmath>

Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM) {
    Matrix& d = r - s;

    Matrix& a = (-GM) * ((d / pow(norm(d), 3)) + (s / pow(norm(s), 3)));

    return a;
}