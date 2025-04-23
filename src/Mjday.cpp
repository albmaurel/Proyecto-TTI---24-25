#include "..\include\Mjday.hpp"

#include <cmath>

double Mjday(int year, int month, int day, int hour=0.0, int minute=0.0, double second=0.0) {
    double jd = 367.0 * year
              - std::floor(7.0 * (year + std::floor((month + 9.0) / 12.0)) * 0.25)
              + std::floor(275.0 * month / 9.0)
              + day + 1721013.5
              + (((second / 60.0 + minute) / 60.0 + hour) / 24.0);

    double Mjd = jd - 2400000.5;

    return Mjd;
}