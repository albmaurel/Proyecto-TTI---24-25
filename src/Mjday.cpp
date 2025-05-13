// $Source$
//----------------------------------------------------------------------
// Mjday
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file Mjday.cpp
 * @brief Calculates Modified Julian Date from calendar date and time
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\Mjday.hpp"


/**
 * @brief Calculates Modified Julian Date from calendar date and time
 * 
 * @param year Year
 * @param month Month
 * @param day Day
 * @param hour Universal time hour
 * @param minute Universal time minute
 * @param second Universal time second
 * 
 * @return Modified Julian Date
 */
double Mjday(int year, int month, int day, int hour=0.0, int minute=0.0, double second=0.0) {
    double jd = 367.0 * year
              - floor(7.0 * (year + floor((month + 9.0) / 12.0)) * 0.25)
              + floor(275.0 * month / 9.0)
              + day + 1721013.5
              + (((second / 60.0 + minute) / 60.0 + hour) / 24.0);

    double Mjd = jd - 2400000.5;

    return Mjd;
}