#include "..\include\GHAMatrix.hpp"

Matrix& GHAMatrix(double Mjd_UT1) {
    
    Matrix& GHAmat = R_z( gast(Mjd_UT1) );

    return GHAmat;
}