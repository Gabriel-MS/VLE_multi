#ifndef RENORMALIZATION_H_INCLUDED
#define RENORMALIZATION_H_INCLUDED

#include "MSA.h"

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include "MSA.h"
#include "numerical.h"

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
            rho = rho/b;
            T = T/b/R*a;
            f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b)+rho*R*T*(log(rho)-1);

            //DIMENSIONLESS!!!************************************************************
            f = b*b/a*f;

            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
            break;

        case 3: //CPA
            rho = rho/b;
            T = T/b/R*a;

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
            //cout << "f = " << f_CPA1 << " / " << f_CPA << " / " << f << endl;
            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);

            //DIMENSIONLESS!!!************************************************************
            f = b*b/a*f;
            break;

        case 4: //MSA
            rho = rho/b;
            T = T/R*a;

            beta = 1/(R*T);
            f = ferg(beta, 0, rho, T, R);

            f = f/a*b;
            break;
    }

   return f;
}

double f0a_function(int EdE, double rho)
{
    double f0a;

    switch(EdE)
    {
        case 1: //SRK
        f0a = -0.5*rho*rho;
        break;


        case 3://CPA
        f0a = -0.5*rho*rho;
        break;
    }

    return f0a;
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

double helmholtz_recursion_long(int EdE, long double f, long double rho, long double a, long double b, double R, double T)
{
    long double fr, f0a;

    switch(EdE)
    {
        case 1: //SRK
            //fr = f + 0.5*a*rho*rho;
            f0a = f0a_function(EdE,rho);

            //DIMENSIONLESS!!!************************************************************
            fr = f - f0a;

            break;

        case 3: //CPA
            //fr = f + 0.5*a*rho*rho;
            f0a = f0a_function(EdE,rho);

            //DIMENSIONLESS!!!************************************************************
            fr = f - f0a;
            break;

        case 4: //MSA
            rho = rho/b;
            T = T/R*a;
            f = f*a/b;

            f0a = f0a_msa(T, 0, rho, 1, 1);
            fr = f - f0a;

            fr = fr/a*b;
            break;
    }

   return fr;
}

double helmholtz_recursion_short(int EdE, long double f, long double rho, long double a, long double b,
                                 int n, double L, long double phi, int sr_type, double R, double T)
{
    long double fr, n2, n2L, c, n2SRK, n2CPA, f0a;

    n2SRK = pow(2,n);
    n2CPA = pow(2,2*n+1);

    switch(EdE)
    {
        case 1: //SRK
        f0a = f0a_function(EdE,rho);

            switch(sr_type)
            {
                case 1: fr = f - f0a*phi/(pow(2,n)); break;
                case 2: fr = f - f0a*phi/(pow(2,2*n+1)); break;
                case 3: fr = f - f0a*phi/(pow(2,2*n-1)); break;

                case 4: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L,2)); break;
                case 5: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L/10,2)); break;
                case 6: fr = f + 0.5*pow(2,-2*n)*rho*rho;
            }
            break;

        case 3: //CPA
            f0a = f0a_function(EdE,rho);

            switch(sr_type)
            {
                case 1: fr = f - f0a*phi/(pow(2,n)); break;
                case 2: fr = f - f0a*phi/(pow(2,2*n+1)); break;
                case 3: fr = f - f0a*phi/(pow(2,2*n-1)); break;

                case 4: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L,2)); break;
                case 5: fr = f + 0.5*phi*a*rho*rho/((pow(2,2*n+1))*pow(L/10,2)); break;
            }
            break;

        case 4: //MSA
            rho = rho/b;
            T = T/R*a;
            f = f*a/b;

            f0a = f0a_msa(T, 0, rho, 2, n);
            fr = f - f0a;

            fr = fr/a*b;
            break;
    }

   return fr;
}

double df_calculation(int w, int n, int Iteration, double width, double Kn, vector<double> rho_vec, VectorXd flv, VectorXd fsv)
{
    std::vector<double> Glv2(1000), Gsv2(1000);
    VectorXd Glv(1000), Gsv(1000), argl(1000), args(1000);
    double suml, sums, aminl, amins, Inl, Ins, rho, al, as, delta_f;
    int t;

                suml = 0;
                sums = 0;

                //Iteração 3 - regra do trapézio para cálculo de I
                t=0;
                aminl = 0;
                amins = 0;

                for(t=0; t<min(w+1,n-w); t++)
                {
                Glv(t) = (flv(w+t) - 2*flv(w) + flv(w-t))/2;
                Gsv(t) = (fsv(w+t) - 2*fsv(w) + fsv(w-t))/2;

                Glv2[t] = exp(-Glv(t)/Kn);
                Gsv2[t] = exp(-Gsv(t)/Kn);

                argl(t) = Glv(t)/Kn;
                args(t) = Gsv(t)/Kn;

                if(argl(t) < aminl) aminl = argl(t);
                if(args(t) < amins) aminl = args(t);
                }

                if(w<(n/2))
                {
                suml = 0.5*(exp(-argl(0)+aminl)+exp(-argl(w)+aminl)); //era 0.25
                sums = 0.5*(exp(-args(0)+amins)+exp(-args(w)+amins));
                }

                else
                {
                suml = 0.5*(exp(-argl(0)+aminl)+exp(-argl(n-w)+aminl));
                sums = 0.5*(exp(-args(0)+amins)+exp(-args(n-w)+amins));
                }

                for(t=1; t<min(w,n-w); t++)
                {
                al = argl(t) - aminl;
                as = args(t) - amins;
                if(al<30) suml = suml + 2.0*exp(-al);
                if(as<30) sums = sums + 2.0*exp(-as);
                }

            Inl = log(width*suml)-aminl;
            Ins = log(width*sums)-amins;


            if(Iteration==2)
            {
                if(w<n/2)
                {
                    Inl = trapezoidal_interval(rho_vec, Glv2, 0, w);
                    Ins = trapezoidal_interval(rho_vec, Gsv2, 0, w);
                }

                else
                {
                    Inl = trapezoidal_interval(rho_vec, Glv2, 0, n-w);
                    Ins = trapezoidal_interval(rho_vec, Gsv2, 0, n-w);
                }

                //Inl = trapezoidal_rule(rho_vec, Glv2);
                //Ins = trapezoidal_rule(rho_vec, Gsv2);
                Inl = log(Inl);
                Ins = log(Ins);
            }

            if(Iteration==3)
            {
                if(w<n/2)
                {
                    Inl = simpson_interval(rho_vec, Glv2, 0, w);
                    Ins = simpson_interval(rho_vec, Gsv2, 0, w);
                }

                else
                {
                    Inl = simpson_interval(rho_vec, Glv2, 0, n-w);
                    Ins = simpson_interval(rho_vec, Gsv2, 0, n-w);
                }

                //Inl = trapezoidal_rule(rho_vec, Glv2);
                //Ins = trapezoidal_rule(rho_vec, Gsv2);
                Inl = log(Inl);
                Ins = log(Ins);
            }

            delta_f = -Kn*(Ins-Inl);

            if(isnan(delta_f) == 1 || isinf(delta_f) == 1) delta_f = 1e-15;

            return delta_f;
}

double critical_exponents()
{
    double Tc, Pc, rho1c, rho2c, rhoc, uc, beta, r2beta, delta, r2delta;
    cout.precision(10);
    std::vector<double> drhoc(1000), rho_d(1000), u_d(1000);
    std::vector<double> exponent(4), T(10), rho1(10), rho2(10), P(10), u(10), lntc(7), lnrhoc(7), lnuc(7), lnpc(7);
    int stop, w, k, r;
    std::string fline;
    int number_lines = 0;
    w = 0;
    k = 0;
    r = 0;

    //Read env data for beta
    double data_beta[12][5];
    ifstream file("../Planilhas de análise/env_exponent.csv");
    cout << "Reading data to calculate critical exponent beta..." << endl;
    for(int row = 0; row < 12; ++row)
    {
    string line;
    getline(file, line);
    if ( !file.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < 5; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> data_beta[row][col];
    }
    }
    file.close();

    //Read env data for delta
    ifstream rfile("../Planilhas de análise/Renormalization.csv");
    std::string rfline;
    int nlines = 0;
    while (std::getline(rfile, rfline))
    {
        ++nlines;
    }

    rfile.close();

    file.open("../Planilhas de análise/Renormalization.csv");
    double data_delta[nlines][8];
    cout << "Reading data to calculate critical exponent delta..." << endl;
    for(int row = 0; row < nlines; ++row)
    {
    string line;
    getline(file, line);
    if ( !file.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < 8; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> data_delta[row][col];
    }
    }
    file.close();

    //Find critical point to calculate beta
    for(int i=0; i<10; i++)
    {
    //i+1 because first line is header
    T[i] = data_beta[i+1][0];
    rho1[i] = data_beta[i+1][1];
    rho2[i] = data_beta[i+1][2];
    P[i] = data_beta[i+1][3];
    u[i] = data_beta[i+1][4];
    }

    Tc = data_beta[11][0];
    rho1c = data_beta[11][1];
    rho2c = data_beta[11][2];
    rhoc = (rho1c+rho2c)/2;
    Pc = data_beta[11][3];
    uc = data_beta[11][4];

    cout << "Difference between rho1c and rho2c = " << fabs(rho1c-rho2c) << " / rhoc = " << rhoc << endl;
    cout << "Tc = " << Tc << endl;
    cout << "Pc = " << Pc << endl;
    cout << "rhoc = " << rhoc << endl;

    //Calculate beta exponent

        //1. Between 0.2% and 1.5% below Critical Temperature

        //2. Between 0.95 and 0.99 Tr

        //3. 1K below Tc to Tc
        for(int i=0; i<8; i++)
        {
        lntc[i] = log(fabs(T[i]-Tc)/Tc);
        lnrhoc[i] = log((rho2[i]-rho1[i])/rhoc);
        }

    exponent = linear_regression(lntc, lnrhoc);
    beta = exponent[0];
    r2beta = exponent[2];

    //Calculate delta exponent
    /*
        //Find data between 0.1 and 0.4 |delta_rho/rhoc|
        int ibegin, iend;
        ibegin = nlines-999;
        iend = nlines;
        for(int i=nlines-999; i<nlines; i++)
        {
        rho_d[i-nlines+999] = data_delta[i][0];
        u_d[i-nlines+999] = data_delta[i][3];
        drhoc[i-nlines+999] = abs(rho_d[i-nlines+999]-rhoc)/rhoc;
        cout << " i / rho / u = " << i-nlines+999 << " / " << rho_d[i-nlines+999] << " / "
             << u_d[i-nlines+999] << " / " << drhoc[i-nlines+999] << endl;
        }

        //Select data of interest
        int ln_size = 0;
        for(int i=0; i<999; i++)
        {
        if(drhoc[i]>=0.1 && drhoc[i]<=0.3) ++ln_size;
        }

        std::vector<double> lndrhoc, lnduc;

        for(int i=0; i<1000; i++)
        {
        if(drhoc[i]>=0.1 && drhoc[i]<=0.3) lndrhoc[i] = log(abs(drhoc[i]/rhoc));
        if(drhoc[i]>=0.1 && drhoc[i]<=0.3) lnduc[i] = log(abs(u_d[i]-uc/uc));
        cout << "drhoc / ln = " << drhoc[i] << " / " << lndrhoc[i] << endl;
        }

        exponent = linear_regression(lndrhoc, lnduc);
        delta = exponent[0];
        r2delta = exponent[2];

        cout << "delta = " << delta << " / " << r2delta << " / " << beta << " / " << r2beta;
    */


    return beta;

}

vector<double> estimate_L_phi(int k, double T)
{
    //variables


    //read experimental data


    //compare actual T with



}

#endif // RENORMALIZATION_H_INCLUDED
