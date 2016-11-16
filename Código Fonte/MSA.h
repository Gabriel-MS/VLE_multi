#ifndef MSA_H_INCLUDED
#define MSA_H_INCLUDED

#include <cmath>

MSA_R(double T, double sigma, double eps, double kB)
{
    double c1, c2, c3, R, Tr;

    Tr = kB*T/eps*6.023e23;

    c1 = 1.1287;
    c2 = -0.05536;
    c3 = 0.0007278;

    R = pow(2,1/6) * pow((1+ pow((1+ (Tr+c2*pow(Tr,2)+c3*pow(Tr,4))/c1),1/2)),-1/6) * sigma;
    return R;
}

MSA_S(double t, double eta)
{
    double S;
    S = pow(1-eta,2)*pow(t,3) + 6*eta*(1-eta)*pow(t,2) + 18*pow(eta,2)*t - 12*eta*(1+2*eta);
    return S;
}

MSA_L(double t, double eta)
{
    double L;
    L = (1+eta/2)*t + 1 + 2*eta;
    return L;
}

MSA_Q(double t, double eta)
{
    double Q, L, S;

    S = MSA_S(t,eta);
    L = MSA_L(t,eta);

    Q = (S+12*eta*L*exp(-t))/pow(1-eta,2)/pow(t,3);
    return Q;
}


helmholtz_MSA(double rho, double T, double sigma, double eps, double kB)
{
    double pi, eta, R, a0, a1, a2, a, beta;
    double a11, a12, a13, g0R, Lz1, Lz2, Qz1, Qz2, z1, z2, a21, a22;
    double k0, k1, k2;

    pi = 3.14159265359;
    k0 = 2.1714*sigma;
    z1 = 2.9637/sigma;
    z2 = 14.0167/sigma;

    R = MSA_R(T, sigma, eps, kB);

    eta = pi*rho*pow(R,3)/6;
    beta = 1/kB/T;

    k1 = k0 * exp(z1*(sigma-R));
    k2 = k0 * exp(z2*(sigma-R));

    g0R = (1+eta/2)/pow(1-eta,2);

    Lz1 = MSA_L(z1*R,eta);
    Lz2 = MSA_L(z2*R,eta);
    Qz1 = MSA_Q(z1*R,eta);
    Qz2 = MSA_Q(z2*R,eta);

    a0 = (4*eta-3*pow(eta,2))/pow(1-eta,2);

    a11 = -12*eta*beta*eps/pow(R,3) * (k1*(Lz1/pow(z1,2)/pow(1-eta,2)/Qz1 - (1+z1*R)/pow(z1,2)) -k2*(Lz2/pow(z2,2)/pow(1-eta,2)/Qz2 - (1+z2*R)/pow(z2,2)));
    a12 = +48*eta*beta*eps * (1/9*pow(sigma/R,6) - 1/3*pow(sigma/R,6));
    a13 = -48*eta*beta*eps*g0R * (1/9*pow(sigma/R,12) - 1/3*pow(sigma/R,6) + 2/9*pow(sigma/R,3));
    a1 = a11 + a12 + a13;

    a21 = -6*eta*pow(beta,2)*pow(eps,2)/pow(R,3) * (pow(k1,2)/2/z1/pow(Qz1,4) + pow(k2,2)/2/z2/pow(Qz2,4) - 2*k1*k2/(z1+z2)/pow(Qz1,2)/pow(Qz2,2));
    a22 = -24*eta*pow(beta,2)*pow(eps,2) * (k1/R/pow(Qz1,2) - k2/R/pow(Qz2,2)) * (1/9*pow(sigma/R,12) - 1/3*pow(sigma/R,6) + 2/9*pow(sigma/R,3));
    a2 = a21 + a22;

    a = a0 + a1 + a2;

    return a;
}


#endif // MSA_H_INCLUDED
