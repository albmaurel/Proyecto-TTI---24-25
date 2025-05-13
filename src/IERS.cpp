// $Source$
//----------------------------------------------------------------------
// IERS
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file IERS.cpp
 * @brief Management of IERS time and polar motion data
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\IERS.hpp"

/**
 * @brief IERS: Management of IERS time and polar motion data
 * 
 * @param eop Matrix containing Earth Orientation Parameters (EOP) data
 * @param Mjd_UTC Modified Julian Date in UTC
 * @param interp Interpolation method ('l' for linear, 'n' for nearest)
 * @return A tuple containing the following values:
 *          - x_pole: x coordinate of the pole in radians    
 *          - y_pole: y coordinate of the pole in radians
 *          - UT1_UTC: UT1-UTC time difference in seconds
 *          - lod: Length of day in seconds
 *          - dpsi: dpsi value in radians
 *          - deps: deps value in radians
 *          - dx_pole: dx coordinate of the pole in radians
 *          - dy_pole: dy coordinate of the pole in radians
 *          - TAI_UTC: TAI-UTC time difference in seconds
 */

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp) {
    

    //return values
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;


    if (interp =='l') {
        // linear interpolation
        double mjd = (floor(Mjd_UTC));
        int i = 0;
        Matrix aux = extract_row(eop,4);
        for (int j = 1; j <= aux.n_column; j++) {
            if (aux(j) == mjd) {
                i = j;
                break;
            }
        }

        Matrix preeop = extract_column(eop,i);

        Matrix nexteop = extract_column(eop,i+1);
        
        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
        y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
        UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
        LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
        dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
        deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
        dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
        dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
        TAI_UTC = preeop(13);
        
        x_pole  = x_pole/SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole  = y_pole/SAT_Const::Arcs;  // Pole coordinate [rad]
        dpsi    = dpsi/SAT_Const::Arcs;
        deps    = deps/SAT_Const::Arcs;
        dx_pole = dx_pole/SAT_Const::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole/SAT_Const::Arcs; // Pole coordinate [rad]
    } else if (interp =='n') {
        double mjd = (floor(Mjd_UTC));
        int i = 0;
        Matrix aux = extract_row(eop,4);
        for (int j = 1; j <= aux.n_column; j++) {
            if (aux(j) == mjd) {
                i = j;
                break;
            }
        }

        Matrix neweop = extract_column(eop,i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = neweop(5)/SAT_Const::Arcs;  // Pole coordinate [rad]
        y_pole  = neweop(6)/SAT_Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = neweop(7);             // UT1-UTC time difference [s]
        LOD     = neweop(8);             // Length of day [s]
        dpsi    = neweop(9)/SAT_Const::Arcs;
        deps    = neweop(10)/SAT_Const::Arcs;
        dx_pole = neweop(11)/SAT_Const::Arcs; // Pole coordinate [rad]
        dy_pole = neweop(12)/SAT_Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = neweop(13);            // TAI-UTC time difference [s]
    }

    return make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}

tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC) {
    return IERS(eop, Mjd_UTC, 'n');
}