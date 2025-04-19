#ifndef __AZELPA__
#define _AZELPA_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"

#include <cmath>

void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif