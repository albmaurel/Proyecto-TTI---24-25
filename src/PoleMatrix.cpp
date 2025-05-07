#include "..\include\PoleMatrix.hpp"

Matrix& PoleMatrix(double xp, double yp) {
    Matrix& polemat = R_y(-xp) * R_x(-yp);

    return polemat;
}