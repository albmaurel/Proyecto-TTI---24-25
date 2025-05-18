#include "..\include\angl.hpp"

double angl(Matrix vec1, Matrix vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;
    double magv1 = norm(vec1);
    double magv2 = norm(vec2);
    double theta  = 0.0;
	double temp;
    if (magv1*magv2 > small^2) {
        temp= dot(vec1,vec2) / (magv1*magv2);
        if (abs( temp ) > 1.0) {
            double sign = temp > 0 ? 1 : -1;
            temp= sign* 1.0;
        }
         theta= acos( temp );
    } else {
        theta= undefined;
    }

    return theta;
}
