// $Source$
//----------------------------------------------------------------------
// gast
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file gast.cpp
 * @brief implementation of the methods auxparam,eop19620101,GGM03S,DE430Coeff
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\global.hpp"


Param AuxParam;

/**
 * @brief Initializes auxiliary parameters for the program.
 *
 * This function sets the values for various auxiliary parameters used in the program.
 * It includes the Modified Julian Date (Mjd_UTC), number of coefficients (n, m),
 * and flags for sun, moon, and planets.
 */
void auxparam() {
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;
    AuxParam.Mjd_TT  = 49746.1170623147;
}

Matrix eopdata;

/**
 * @brief Reads Earth Orientation Parameters (EOP) data from a file and stores it in a matrix.
 *
 * @param c The number of columns (data entries) to read from the file.
 */
void eop19620101(int c) {
    eopdata = zeros(13, c);

    FILE *fid = fopen("..\\data\\eop19620101.txt", "r");
    if (fid== NULL) {
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 1; i <= c; i++) {
        fscanf(fid,"%lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &(eopdata(1,i)),&(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),&(eopdata(5,i)),
            &(eopdata(6,i)),&(eopdata(7,i)),&(eopdata(8,i)),&(eopdata(9,i)),&(eopdata(10,i)),
            &(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
    }
    fclose(fid);
}

Matrix Cnm;
Matrix Snm;

/**
 * @brief Reads GGM03S gravitational coefficients from a file and stores them in matrices.
 */
void GGM03S() {
    
    Cnm = zeros(181, 181);
    Snm = zeros(181, 181);
    
    FILE *fid = fopen("..\\data\\GGM03S.txt", "r");
    if (fid== NULL) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for(int n = 0; n <= 180; n++) {
        for(int m = 0; m <= n; m++) {
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",
                &aux,&aux,&(Cnm(n+1, m+1)),&(Snm(n+1, m+1)),
                &aux,&aux);
        }
    }
    fclose(fid);
}

Matrix PC;

/**
 * @brief Reads DE430 planetary coefficients from a file and stores them in a matrix.
 */
void DE430Coeff() {
    PC = zeros(2285, 1020);
	FILE *fid = fopen("../data/DE430Coeff.txt","r");

	if(fid== NULL) {
		printf("Fail open DE430Coeff.txt file\n");
		exit(EXIT_FAILURE);
	}
	double aux;
	for (int n = 1; n <= 2285; n++) {
		for(int m=1;m<=1020;m++){
				fscanf(fid, "%lf ",&(PC(n, m)));
			}
		}
	fclose(fid);
}

