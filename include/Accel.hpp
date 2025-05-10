#ifndef _ACCEL_
#define _ACCEL_

#include "../include/matrix.hpp"
#include "../include/SAT_Const.hpp"
#include "../include/global.hpp"
#include "../include/AccelPointMass.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/JPL_Eph_DE430.hpp"
#include "../include/Mjday_TDB.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/timediff.hpp"
#include "../include/IERS.hpp"

Matrix& Accel(double x, Matrix Y);

#endif