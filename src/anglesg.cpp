// $Source$
//----------------------------------------------------------------------
// anglesg
//----------------------------------------------------------------------
// Proyecto TT1
//
/** @file anglesg.cpp
 * @brief This file contains the implementation of the anglesg function.
 * 
 * @author Alberto Maurel Mendizábal
 * @bug No known bugs.
 */
//----------------------------------------------------------------------
#include "..\include\anglesg.hpp"

/**
 * @brief Computes the geocentric inertial position and velocity vectors using the Gauss angles-only method.
 *
 * This function implements the Gauss method for initial orbit determination using three optical observations
 * (azimuth and elevation angles) at three different times, along with the observer's position vectors at those times.
 * It transforms the observed directions to the inertial frame, solves for the range to the object, and computes
 * the position and velocity vectors in the geocentric inertial frame.
 *
 * @param az1 Azimuth angle at the first observation (radians).
 * @param az2 Azimuth angle at the second observation (radians).
 * @param az3 Azimuth angle at the third observation (radians).
 * @param el1 Elevation angle at the first observation (radians).
 * @param el2 Elevation angle at the second observation (radians).
 * @param el3 Elevation angle at the third observation (radians).
 * @param Mjd1 Modified Julian Date of the first observation (UTC).
 * @param Mjd2 Modified Julian Date of the second observation (UTC).
 * @param Mjd3 Modified Julian Date of the third observation (UTC).
 * @param Rs1 Observer's position vector at the first observation (ECEF frame).
 * @param Rs2 Observer's position vector at the second observation (ECEF frame).
 * @param Rs3 Observer's position vector at the third observation (ECEF frame).
 * @return A tuple containing:
 *         - r2: The geocentric inertial position vector at the second observation (Matrix&).
 *         - v2: The geocentric inertial velocity vector at the second observation (Matrix&).
 */
tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3) {
	//Las variables que están inicializadas aquí son aquellas que entran por primera vez en un bucle o en un método.
    Matrix L1(3); Matrix L2(3); Matrix L3(3);
    Matrix Lm1(3), Lm2(3), Lm3(3);
	Matrix z(3);
	double x[10], y[10];
    Matrix r1(3), r2(3), r3(3);
	
    double magr1, magr2, magr3;
    Matrix v2(3);
    double theta, theta1, copa;
    string error;
    double p, a, e, i, Omega, omega, M;
    double f1, f3, g1, g3;
    double rdot, udot;
    double H1, H2, H3;
    double G1, G2, G3;
    double D1, D2, D3;
    double tausqr;
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;

    L1(1)=cos(el1)*sin(az1); L1(2)= cos(el1)*cos(az1); L1(3)= sin(el1);;
    L2(1)=cos(el2)*sin(az2); L2(2)= cos(el2)*cos(az2); L2(3)= sin(el2);;
    L3(1)=cos(el3)*sin(az3); L3(2)= cos(el3)*cos(az3); L3(3)= sin(el3);;

    auto [lon1, lat1, h1] = Geodetic(Rs1);
    auto [lon2, lat2, h2] = Geodetic(Rs2);
    auto [lon3, lat3, h3] = Geodetic(Rs3);

    Matrix M1 = LTC(lon1, lat1);
    Matrix M2 = LTC(lon2, lat2);
    Matrix M3 = LTC(lon3, lat3);    

    // body-fixed system
    Matrix Lb1 = transpose(M1)*transpose(L1);
    Matrix Lb2 = transpose(M1)*transpose(L2);
    Matrix Lb3 = transpose(M1)*transpose(L3);;

    // mean of date system (J2000)
    
    double Mjd_UTC = Mjd1;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');	timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P = PrecMatrix(SAT_Const::MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Lm1 = transpose(E)*Lb1;
    Rs1 = transpose(E)*Rs1;
    
    Mjd_UTC = Mjd2;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
	timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(SAT_Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Lm2 = transpose(E)*Lb2;
    Rs2 = transpose(E)*Rs2;
    
    Mjd_UTC = Mjd3;
    tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
	timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(SAT_Const::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Lm3 = transpose(E)*Lb3;
    Rs3 = transpose(E)*Rs3;

    // geocentric inertial position
    double tau1 = (Mjd1-Mjd2)*86400;
    double tau3 = (Mjd3-Mjd2)*86400;

    double a1 = tau3/(tau3-tau1);
    double a3 =-tau1/(tau3-tau1);

    double b1 = tau3/(6*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau3,2.0));
    double b3 =-tau1/(6*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau1,2.0));

    Matrix f(3,3);
    assign_column(f,transpose(Lm1),1);
    assign_column(f,transpose(Lm2),2);
    assign_column(f,transpose(Lm3), 3);
    Matrix  b(3,3);
    assign_column(b,transpose(Rs1),1);
    assign_column(b,transpose(Rs2),2);
    assign_column(b,transpose(Rs2),3);

    Matrix D = inv(f)*b;

    double d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    double d2s = D(2,1)*b1+D(2,3)*b3;

    double Ccye = 2*dot(Lm2,Rs2);
    
    double poly[10] = {
        1.0,  // R2^8... polynomial
        0.0,
        -(pow(d1s,2) + d1s*Ccye + pow(norm(Rs2),2)),
        0.0,
        0.0,
        -SAT_Const::GM_Earth*(d2s*Ccye + 2*d1s*d2s),
        0.0,
        0.0,
        -pow(SAT_Const::GM_Earth,2)*pow(d2s,2)
    };
    real_poly_roots(poly, 8, x, y); // es lo que se modifica debido al cambio en la función
    double bigr2= -99999990.0;

    for (int j=1; j <= 8; j++) {
        if (( x[j] > bigr2 ) && ( y[j]==0))  {
            bigr2= x[j];
        }
    }

    double u = SAT_Const::GM_Earth/(pow(bigr2,3));
    
    double C1 = a1+b1*u;
    double C2 = -1;
    double C3 = a3+b3*u;

    z(1)=C1;
    z(2)=C2;
    z(3)=C3;

    Matrix temp = -D*transpose(z);
    double rho1 = temp(1)/(a1+b1*u);
    double rho2 = -temp(2);
    double rho3 = temp(3)/(a3+b3*u);

    double rhoold1 = rho1;
    double rhoold2 = rho2;
    double rhoold3 = rho3;

    rho2 = 99999999.9;
    int ll   = 0;

    while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
        ll = ll + 1;
        rho2 = rhoold2;
        
        r1 = Rs1+Lm1*rho1;
        r2 = Rs2+Lm2*rho2;
        r3 = Rs3+Lm3*rho3;
        
        magr1 = norm(r1);
        magr2 = norm(r2);
        magr3 = norm(r3);
        
        tie(v2, theta,theta1,copa,error) = gibbs(r1,r2,r3);
        
        if ( (error!="          ok") & (copa < M_PI/180) ) {
            tie(v2,theta,theta1,copa,error) = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);
        }
		//[r2,v2]
        Matrix temp4 = union_vector(transpose(r2), transpose(v2));
        auto [p, a, e, i, Omega, omega, M] = elements (temp4);
		
        
        if ( ll <= 8)  {
            u = SAT_Const::GM_Earth/pow(magr2,3);
            rdot= dot(r2,v2)/magr2;
            udot= (-3*SAT_Const::GM_Earth*rdot)/(pow(magr2,4));
            
            tausqr= tau1*tau1;
            f1=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau1 
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau1;
            g1= tau1 - (1/6)*u*tau1*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau1 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau3*tau3;
            f3=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau3 
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau3;
            g3= tau3 - (1/6)*u*tau3*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau3 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
        } else {
            
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );
            
            f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
            f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
        }
        
        C1 = g3/(f1*g3-f3*g1);
        C2 = -1;
        C3 =-g1/(f1*g3-f3*g1);
        
        H1 = SAT_Const::GM_Earth*tau3/12;
        H3 =-SAT_Const::GM_Earth*tau1/12;
        H2 = H1-H3;
        
        G1 = -tau3/(tau1*(tau3-tau1));
        G3 = -tau1/(tau3*(tau3-tau1));
        G2 = G1-G3;
        
        D1 = G1+H1/pow(magr1,3);
        D2 = G2+H2/pow(magr2,3);
        D3 = G3+H3/pow(magr3,3);
		//temp = -[D1 D2 D3]*[C1 C2 C3]';
        Matrix temp1(3), temp3(3);
        temp1(1)=D1; temp1(2)=D2; temp1(3)=D3;
        temp3(1)=C1; temp3(2)=C2; temp3(3)=C3;
        double temp2 = (-temp1*transpose(temp3))(1,1);
        rhoold1 = temp2/(a1+b1*u);
        rhoold2 = -temp2;
        rhoold3 = temp2/(a3+b3*u);
        
        r1 = Rs1+Lm1*rhoold1;
        r2 = Rs2+Lm2*rhoold2;
        r3 = Rs3+Lm3*rhoold3;
        
    }

    r1 = Rs1+Lm1*rho1;
    r2 = Rs2+Lm2*rho2;
    r3 = Rs3+Lm3*rho3;

    return tie(r2, v2);
}