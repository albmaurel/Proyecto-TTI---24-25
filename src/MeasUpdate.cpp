// $Source$
//----------------------------------------------------------------------
// MeasUpdate
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file MeasUpdate.cpp
 * @brief Performs a measurement update step in a Kalman filter.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\MeasUpdate.hpp"

/**
 * @brief Performs a measurement update step in a Kalman filter.
 * 
 * @param x The current state vector.
 * @param z The measurement value.
 * @param g The predicted measurement value.
 * @param s The standard deviation of the measurement noise.
 * @param G The measurement matrix.
 * @param P The current state covariance matrix.
 * @param n The dimension of the state vector.
 * @return A tuple containing:
 *         - The Kalman gain matrix (Matrix& K).
 *         - The updated state vector (Matrix& x2).
 *         - The updated state covariance matrix (Matrix& P2).
 */
tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n) {
    double m = 1;
    Matrix Inv_W = zeros(m,m);

    Inv_W(1,1) = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    Matrix& K = P*transpose(G)*inv(Inv_W+G*P*transpose(G));

    // State update
    Matrix& x2 = x + K*(z-g);

    //Covariance update
    Matrix& P2 = (eye(n)-K*G)*P;

    return tie(K, x2, P2);
}
