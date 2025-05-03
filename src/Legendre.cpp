#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <iostream>

#include <cmath>
#include <tuple>

tuple<Matrix, Matrix> Legendre(int n, int m, double fi)
{
    Matrix pnm = zeros(n+1, m+1);
    Matrix dpnm = zeros(n+1, m+1);

    pnm(1, 1) = 1;
    dpnm(1, 1) = 0;
    pnm(2,2)=sqrt(3.0)*cosf(fi);
    dpnm(2,2)=-sqrt(3.0)*sinf(fi);
    //diagonal coefficients
    for(size_t i = 2; i<= n; i++){
        pnm(i+1,i+1)= sqrt((2.0*i+1)/(2*i))*cosf(fi)*pnm(i,i);
        dpnm(i+1,i+1)= sqrt((2.0*i+1)/(2*i))*((cosf(fi)*dpnm(i,i))-(sinf(fi)*pnm(i,i)));
    }
    // horizontal first step coefficients
    for(size_t i = 1; i<= n; i++){ 
        pnm(i+1,i)= sqrt(2.0*i+1)*sinf(fi)*pnm(i,i);
    
        dpnm(i+1,i)= sqrt(2.0*i+1)*((cosf(fi)*pnm(i,i))+(sinf(fi)*dpnm(i,i)));
    }
    // horizontal second step coefficients
    int j=0;
    int k=2;
    while(1){
        for(size_t i = k; i<= n; i++){         
            pnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sinf(fi)*pnm(i,j+1))
            -(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-1,j+1)));
            dpnm(i+1,j+1)=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sinf(fi)*dpnm(i,j+1))
                +(sqrt(2*i-1)*cosf(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2.0*i-3))*dpnm(i-1,j+1)));
        }
        j = j+1;
        k = k+1;
        if (j>m){
            break;
        }
    }
    return make_tuple(pnm, dpnm);
}