#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"

#include <cmath>
#include <tuple>
tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp);

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC);

#endif