#include "..\include\matrix.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\SAT_Const.hpp" 


double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT - SAT_Const::MJD_J2000) / 36525.0;

    double arcsec = 84381.448 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T;

    return SAT_Const::Rad * (arcsec / 3600.0);  
}