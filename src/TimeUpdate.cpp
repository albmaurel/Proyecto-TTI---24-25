#include "..\include\matrix.hpp"
#include "..\include\TimeUpdate.hpp"


Matrix TimeUpdate(const Matrix& P, const Matrix& Phi, double Qdt) {
	Matrix PhiT = transpose(Phi);
	Matrix updatedP = Phi * P * PhiT + Qdt;
	
	return updatedP;
}