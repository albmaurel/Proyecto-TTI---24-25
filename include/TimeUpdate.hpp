#ifndef _TIME_UPDATE_
#define _TIME_UPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix TimeUpdate(const Matrix& P, const Matrix& Phi, double Qdt = 0.0);

#endif