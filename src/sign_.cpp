#include "..\include\sign_.hpp"

#include <cmath>

double sign_(double a, double b) {
    return (b >= 0.0) ? abs(a) : -abs(a);
}