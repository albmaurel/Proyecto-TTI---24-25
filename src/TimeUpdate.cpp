#include "..\include\matrix.hpp"
#include "..\include\TimeUpdate.hpp"


Matrix TimeUpdate( Matrix& P,  Matrix& Phi, double Qdt) {
	Matrix PhiT = transpose(Phi);
	Matrix updatedP = Phi * P * PhiT + Qdt;
	
	return updatedP;
}