#include "..\include\matrix.hpp"
#include "..\include\AccelPointMass.hpp"
#include <cmath>

Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM) {
    Matrix& d = r - s;

    Matrix temp = (d / pow(norm(d), 3)) + (s / pow(norm(s), 3));
    Matrix& a = temp * (-GM);

    return a;
}