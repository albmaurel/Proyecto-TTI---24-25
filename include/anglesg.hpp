#ifndef _ANGLESG_
#define _ANGLESG_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\elements.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\global.hpp"
#include "..\include\IERS.hpp"
#include "..\include\geodetic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\timediff.hpp"
#include "..\include\hgibbs.hpp"
#include "..\include\gibbs.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\rpoly.hpp"
#include <tuple>

tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3);

#endif