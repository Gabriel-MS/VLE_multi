#ifndef RENORMALIZATION_H_INCLUDED
#define RENORMALIZATION_H_INCLUDED

#include "MSA.h"

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;

double helmholtz_repulsive(int EdE, double R, double T, long double rho, long double a, long double b, VectorXd X, VectorXd x,
                           double sigma, double eps, double kB)
{
    double f, f_CPA, f1;
    VectorXd f_CPA1(8);
    MatrixXd one_4c(8,2);

    switch(EdE)
    {
        case 1: //SRK
            //f = rho*R*T*(log(rho/(1-rho*b))-1)-rho*a/b*log(1+rho*b);
            f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b)+rho*R*T*(log(rho)-1);
            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
            break;

        case 3: //CPA
            one_4c << 1, 0,
                      1, 0,
                      1, 0,
                      1, 0,
                      0, 1,
                      0, 1,
                      0, 1,
                      0, 1;

            f_CPA1 = X.array().log()-0.5*X.array()+0.5;
            //f_CPA = (one_4c*x).transpose()*f_CPA1;
            f_CPA = x.transpose()*(f_CPA1.transpose()*one_4c);

            f = rho*R*T*(log(rho/(1-rho*b))-1)-rho*a/b*log(1+rho*b)+rho*R*T*f_CPA;
            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
            break;

        case 4: //MSA
            f1 = helmholtz_MSA(rho, T, sigma, eps, kB);
            f = rho*R*T*f1 + rho*R*T*(log(rho)-1);
            break;
    }

   return f;
}


double helmholtz_recursion_long(int EdE, long double f, long double rho, long double a)
{
    long double fr;

    switch(EdE)
    {
        case 1: //SRK
            fr = f + 0.5*a*rho*rho;
            break;

        case 3: //CPA
            fr = f + 0.5*a*rho*rho;
            break;

        case 4: //MSA
            fr = f + a*rho*rho;
    }

   return fr;
}

double helmholtz_recursion_short(int EdE, long double f, long double rho, double a, int n, double L, long double phi, int sr_type)
{
    long double fr, n2, n2L, c, n2SRK, n2CPA;

    n2SRK = pow(2,n);
    n2CPA = pow(2,2*n+1);

    switch(EdE)
    {
        case 1: //SRK

            switch(sr_type)
            {
                case 1: fr = f + 0.5*phi*a*rho*rho/(pow(2,n)); break;
                case 2: fr = f + 0.5*phi*a*rho*rho/(pow(2,2*n+1)); break;
                case 3: fr = f + 0.5*phi*a*rho*rho/(pow(2,2*n-1)); break;
                case 4: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L,2)); break;
                case 5: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L/10,2)); break;
            }
            break;

        case 3: //CPA
            //fr = f + 0.5*a*rho*rho*phi/n2CPA;

            switch(sr_type)
            {
                case 1: fr = f + 0.5*phi*a*rho*rho/(pow(2,n)); break;
                case 2: fr = f + 0.5*phi*a*rho*rho/(pow(2,2*n+1)); break;
                case 3: fr = f + 0.5*phi*a*rho*rho/(pow(2,2*n-1)); break;
                case 4: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L,2)); break;
                case 5: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L/10,2)); break;
            }
            break;

        case 4: //MSA

            switch(sr_type)
            {
                case 1: fr = f + phi*a*rho*rho/(pow(2,n)); break;
                case 2: fr = f + phi*a*rho*rho/(pow(2,2*n+1)); break;
                case 3: fr = f + phi*a*rho*rho/(pow(2,2*n-1)); break;
                case 4: fr = f + phi*a*rho*rho/((pow(2,2*n+1))*pow(L,2)); break;
                case 5: fr = f + phi*a*rho*rho/((pow(2,2*n+1))*pow(L/10,2)); break;
            }
            break;
    }

   return fr;
}


#endif // RENORMALIZATION_H_INCLUDED
