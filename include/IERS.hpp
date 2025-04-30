#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"

#include <cmath>

void IERS( Matrix& eop, double Mjd_UTC, char interp,
          double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

#endif