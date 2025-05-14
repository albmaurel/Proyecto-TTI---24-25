#ifndef _HGIBBS_
#define _HGIBBS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include <cmath>
#include <tuple>

tuple<Matrix&, double, double, double, string> hgibbs(Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3);

#endif