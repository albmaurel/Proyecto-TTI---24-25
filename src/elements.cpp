#include "..\include\elements.hpp"

/**
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits
 * 
 * @param y State  (x,y,z,vx,vy,vz)
 * 
 * @return tuple containing:
 *         - double p: Semilatus rectum [m]
 *         - double a: Semimajor axis
 *         - double e: Eccentricity
 *         - double i: Inclination (radians)
 *         - double Omega: Longitude of the ascending node [rad]
 *         - double omega: Argument of pericenter [rad]
 *         - double M: Mean anomaly [rad]
 */

tuple<double, double, double, double, double, double, double> elements(Matrix y) {
        
    double pi2 = SAT_Const::pi2;
    
    Matrix r =extract_vector(y,1,3);                                        // Position
    Matrix v =extract_vector(y,4,6);                                        // Velocity
    
    Matrix h = cross(r,v);                                    // Areal velocity
    double magh = norm(h);
    double p = magh*magh/SAT_Const::GM_Earth;
    double H = norm(h);
    
    double Omega = atan2 ( h(1), -h(2) );                     // Long. ascend. node 
    Omega = fmod(Omega,pi2);
    if (Omega < 0) {
        Omega += pi2;
    }
    double i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); // Inclination        
    double u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    // Arg. of latitude   
    
    double R  = norm(r);                                      // Distance           
    
    double a = 1.0/(2.0/R-dot(v,v)/SAT_Const::GM_Earth);           // Semi-major axis    
    
    double eCosE = 1.0-R/a;                                      // e*cos(E)           
    double eSinE = dot(r,v)/sqrt(SAT_Const::GM_Earth*a);       // e*sin(E)           
    
    double e2 = eCosE*eCosE +eSinE*eSinE;
    double e  = sqrt(e2);                                      // Eccentricity 
    double E  = atan2(eSinE,eCosE);                            // Eccentric anomaly  
    
    double M  = fmod(E-eSinE,pi2);                             // Mean anomaly
    if (M < 0) {
        M+=pi2;
    }
    
    double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);           // True anomaly
    
    double omega = fmod(u-nu,pi2);                             // Arg. of perihelion 
    if (omega < 0) {
        omega += pi2;
    }
    
    return tie(p, a, e, i, Omega, omega, M);
    
}