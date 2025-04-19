#include "..\include\IERS.hpp"
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"

#include <cmath>

#include <cmath>

void IERS(const Matrix& eop, double Mjd_UTC, char interp,
          double& x_pole, double& y_pole, double& UT1_UTC, double& LOD,
          double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {

    if (interp != 'l' && interp != 'n') {
        interp = 'n'; 
    }

    double mjd = std::floor(Mjd_UTC);
    int i = -1; /

    for (int j = 1; j <= eop.n_column; ++j) {
        if (mjd == eop(4, j)) {
            i = j;
            break;
        }
    }

    if (i == -1) {
        cout << "Error: Mjd_UTC not found in eop data.\n";
		exit(EXIT_FAILURE);
    }

    if (interp == 'l') {
        if (i + 1 > eop.n_column) {
            cout << "Error: Not enough eop data for linear interpolation.\n";
			exit(EXIT_FAILURE);
        }

        double mfme = 1440.0 * (Mjd_UTC - std::floor(Mjd_UTC));
        double fixf = mfme / 1440.0;

        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = eop(5, i) + (eop(5, i + 1) - eop(5, i)) * fixf;
        y_pole = eop(6, i) + (eop(6, i + 1) - eop(6, i)) * fixf;
        UT1_UTC = eop(7, i) + (eop(7, i + 1) - eop(7, i)) * fixf;
        LOD = eop(8, i) + (eop(8, i + 1) - eop(8, i)) * fixf;
        dpsi = eop(9, i) + (eop(9, i + 1) - eop(9, i)) * fixf;
        deps = eop(10, i) + (eop(10, i + 1) - eop(10, i)) * fixf;
        dx_pole = eop(11, i) + (eop(11, i + 1) - eop(11, i)) * fixf;
        dy_pole = eop(12, i) + (eop(12, i + 1) - eop(12, i)) * fixf;
        TAI_UTC = eop(13, i);
        x_pole = x_pole / Arcs;
        y_pole = y_pole / Arcs;
        dpsi = dpsi / Arcs;
        deps = deps / Arcs;
        dx_pole = dx_pole / Arcs;
        dy_pole = dy_pole / Arcs;
    } else if (interp == 'n') {
        // No interpolation
        x_pole = eop(5, i) / Arcs;
        y_pole = eop(6, i) / Arcs;
        UT1_UTC = eop(7, i);
        LOD = eop(8, i);
        dpsi = eop(9, i) / Arcs;
        deps = eop(10, i) / Arcs;
        dx_pole = eop(11, i) / Arcs;
        dy_pole = eop(12, i) / Arcs;
        TAI_UTC = eop(13, i);
    }
}