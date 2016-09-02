#ifndef RENORMALIZATION_H_INCLUDED
#define RENORMALIZATION_H_INCLUDED

double helmholtz_density(int EdE, double R, double T, double rho, double a, double b)
{
    double f;

    switch(EdE)
    {
        case 1: //SRK
            f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
            break;
    }

   return f;
}



double helmholtz_repulsive(int EdE, double R, double T, double rho, double a, double b)
{
    double f;

    switch(EdE)
    {
        case 1: //SRK
            f = rho*R*T*(log(rho/(1-rho*b))-1)-rho*a/b*log(1+rho*b);
            break;
    }

   return f;
}


double helmholtz_recursion_long(int EdE, double f, double rho, double a)
{
    double fr;

    switch(EdE)
    {
        case 1: //SRK
            fr = f + 0.5*a*rho*rho;
            break;
    }

   return fr;
}

double helmholtz_recursion_short(int EdE, double f, double rho, double a, int n, double L)
{
    double fr, n2, n2L, c;

    n2 = pow(2,n);
    n2L = pow(n2*L,3);
    c = 0.5;

    switch(EdE)
    {
        case 1: //SRK
            fr = f + 0.5*c*a*rho*rho/n2; //0.5 is for alkanes
            break;
    }

   return fr;
}


#endif // RENORMALIZATION_H_INCLUDED
