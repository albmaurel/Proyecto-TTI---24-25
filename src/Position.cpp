#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Position.hpp"


#include <cmath>

Matrix Position(double lon, double lat, double h) {
    double R_equ = R_Earth;
    double f     = f_Earth;

    double e2     = f * (2.0 - f);           // Square of eccentricity
    double CosLat = std::cos(lat);			//Co)sine of geodetic latitude
    double SinLat = std::sin(lat);
	
	//Position vector
    double N = R_equ / std::sqrt(1.0 - e2 * SinLat * SinLat);

    Matrix r = zeros(3);

    r(1) = (N + h) * CosLat * std::cos(lon);
    r(2) = (N + h) * CosLat * std::sin(lon);
    r(3) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}