#ifndef RENORMALIZATION_H_INCLUDED
#define RENORMALIZATION_H_INCLUDED

#include "MSA.h"

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include "MSA.h"

using namespace Eigen;

double helmholtz_repulsive(int EdE, double R, double T, long double rho, long double a, long double b, VectorXd X, VectorXd x,
                           double sigma, double eps, double kB)
{
    double f, f_CPA, f1, beta;
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
            f_CPA = (one_4c*x).transpose()*f_CPA1;

            f = rho*R*T*(log(rho/(1-rho*b))-1)-rho*a/b*log(1+rho*b)+rho*R*T*f_CPA;
            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
            break;

        case 4: //MSA
            beta = 1/(R*T);
            f = ferg(beta, 0, rho, T, R);
            break;
    }

   return f;
}

double f0a_msa(double T, double rho1, double rho2, int type, int n)
{
    double f0a;
    double b11, b22, b12, a11, a22, a12, xlam11, xlam12, xlam22, del11, del22, del12;
    double b111, b112, b113, a111, a112, a113, b221, b222, b223, a221, a222, a223, pi, cnst;
    double sig1, sig2, sig12, cnst1, xkmix, xlmix;
    vector<double> data(7);

pi = 3.14159265359796;
cnst = pi/6;
cnst1 = 6/pi;

data = d_msa(1);
a111 = data[0];
a112 = data[1];
a113 = data[2];
b111 = data[3];
b112 = data[4];
b113 = data[5];
xlam11 = data[6];

data = d_msa(2);
a221 = data[0];
a222 = data[1];
a223 = data[2];
b221 = data[3];
b222 = data[4];
b223 = data[5];
xlam22 = data[6];

xlmix = 0.00000E-03;
xkmix = 0.18504E+00;
xlam12 = (xlam11+xlam22)/2;

b11 = b111 + b112/T + b113/pow(T,2);
b22 = b221 + b222/T + b223/pow(T,2);
a11 = (a111+pow((a112/T),a113));
a22 = (a221+pow((a222/T),a223));

sig1 = pow((b11/cnst),(1.0/3.0));
sig2 = pow((b22/cnst),(1.0/3.0));
sig12 = 0.5*(sig1+sig2)*(1.0-xlmix);
b12 = cnst*pow(sig12,3.0);
a12 = pow((a11*a22),0.5)*(1.0-xkmix);

    switch(type)
    {
    case 1: //long
        f0a = -4.0*(a11*b11*pow(xlam11,3)*rho1*rho1 +a22*b22*pow(xlam22,3)*rho2*rho2 +2.0*a12*b12*pow(xlam12,3)*rho1*rho2);
        break;

    case 2: //short
        del11 = 1.16;
        del22 = 2.62;
        del12 = 10.6;
        f0a = pow(2.0,-2*n)*4.0 *(-del11*a11*b11*pow(xlam11,3)*rho1*rho1 - del22*a22*b22*pow(xlam22,3)*rho2*rho2 -2.0*del12*a12*b12*pow(xlam12,3)*rho1*rho2);
        break;
    }

    return f0a;
}

double helmholtz_recursion_long(int EdE, long double f, long double rho, long double a)
{
    long double fr, f0a;

    switch(EdE)
    {
        case 1: //SRK
            fr = f + 0.5*a*rho*rho;
            break;

        case 3: //CPA
            fr = f + 0.5*a*rho*rho;
            break;

        case 4: //MSA
            f0a = f0a_msa(500, 0, rho, 1, 1);
            fr = f - f0a;
            break;
    }

   return fr;
}

double helmholtz_recursion_short(int EdE, long double f, long double rho, double a, int n, double L, long double phi, int sr_type)
{
    long double fr, n2, n2L, c, n2SRK, n2CPA, f0a;

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
            f0a = f0a_msa(500, 0, rho, 2, n);
            fr = f - f0a;
            break;
    }

   return fr;
}


#endif // RENORMALIZATION_H_INCLUDED
