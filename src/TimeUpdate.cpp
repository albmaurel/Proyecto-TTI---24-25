// $Source$
//----------------------------------------------------------------------
// TimeUpdate
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file TimeUpdate.cpp
 * @brief Performs the time update step for a covariance matrix in a Kalman filter.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\TimeUpdate.hpp"


/**
 * @brief Performs the time update step for a covariance matrix in a Kalman filter.
 * 
 * 
 * @param P Reference to the current covariance matrix to be updated.
 * @param Phi Reference to the state transition matrix.
 * @param Qdt Process noise covariance scaled by the time step.
 * @return The updated covariance matrix after the time update step.
 */
Matrix TimeUpdate( Matrix& P,  Matrix& Phi, double Qdt) {
	Matrix PhiT = transpose(Phi);
	Matrix updatedP = Phi * P * PhiT + Qdt;
	
	return updatedP;
}