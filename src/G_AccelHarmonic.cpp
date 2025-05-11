#include "..\include\G_AccelHarmonic.hpp"

Matrix& G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max) {
    double d = 1.0;   // Position increment [m]

    Matrix& G = zeros(3,3);
    Matrix& dr = zeros(3,1);
    Matrix da;
    // Gradient
    for (int i=1; i <= 3; i++) {
        // Set offset in i-th component of the position vector
        dr = zeros(3,1);
        dr(i) = d;
        // Acceleration difference
        da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) -
        AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis

        assign_column(G,transpose((da/d)),i);
    }
    return G;
}
