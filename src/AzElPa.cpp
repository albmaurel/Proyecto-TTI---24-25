#include "..\include\matrix.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\SAT_Const.hpp"


#include <cmath>

void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {
	double pi2 = SAT_Const::pi2;
    double rho = std::sqrt(s(1) * s(1) + s(2) * s(2));

    // Angles
    Az = std::atan2(s(1), s(2)); // Azimuth

    if (Az < 0.0) {
        Az += pi2;
    }

    El = std::atan(s(3) / rho); // Elevation

    // Partials
    dAds(1) = s(2) / (rho * rho);
    dAds(2) = -s(1) / (rho * rho);
    dAds(3) = 0.0;

    dEds(1) = -s(1) * s(3) / rho;
    dEds(2) = -s(2) * s(3) / rho;
    dEds(3) = rho;

    // Normalize dEds by the dot product of s
    double dot_product = dot(s,s);
    dEds /= dot_product;
}