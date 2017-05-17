#ifndef RENORMALIZATION_H_INCLUDED
#define RENORMALIZATION_H_INCLUDED

#include "MSA.h"

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include "MSA.h"
#include "numerical.h"
#include "bicubic.h"
#include "math.h"
#include "EdE.h"

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

        case 5: //SRK for binary, only residual part
            //f = rho*R*T*(log(rho/(1-rho*b))-1)-rho*a/b*log(1+rho*b);
            rho = rho/b;
            T = T/b/R*a;
            f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);

            //DIMENSIONLESS!!!************************************************************
            f = b*b/a*f;

            //f = -rho*R*T*log(1-rho*b)-rho*a/b*log(1+rho*b);
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
    std::vector<double> Glv2(n), Gsv2(n);
    VectorXd Glv(n), Gsv(n), argl(n), args(n);
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

double critical_exponents(int EdE)
{
    double Tc, Pc, rho1c, rho2c, rhoc, uc, r2beta, delta, r2delta, beta2;
    cout.precision(10);
    std::vector<double> drhoc(1000), rho_d(1000), u_d(1000);
    std::vector<double> exponent(3), T(10), rho1(10), rho2(10), P(10), u(10), lntc(7), lnrhoc(7), lnuc(7), lnpc(7);
    int stop, w, k, r;
    std::string fline;
    int number_lines = 0;
    w = 0;
    k = 0;
    r = 0;

    //Read env data for beta
    double data[12][5];
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
    convertor >> data[row][col];
    }
    }
    file.close();
/*
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
*/
    //Find critical point to calculate beta
    for(int i=0; i<10; i++)
    {
    //i+1 because first line is header
    T[i] = data[i+1][0];
    rho1[i] = data[i+1][1];
    rho2[i] = data[i+1][2];
    P[i] = data[i+1][3];
    u[i] = data[i+1][4];
    }

    Tc = data[11][0];
    rho1c = data[11][1];
    rho2c = data[11][2];
    rhoc = (rho1c+rho2c)/2;
    Pc = data[11][3];
    uc = data[11][4];

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
    //beta2 = exponent[0];
    beta2 = linear_regression_angular(lntc, lnrhoc);
    r2beta = exponent[2];

    cout << "beta inside = " << beta2 << endl;
    cout << "r2beta inside = " << r2beta << endl;

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

    cout << "almost out of critical exponents = " << beta2 << endl;
    return beta2;
}

void beta_exponent(double *beta, double Tc, double rhoc)
{
    double Pc, rho1c, rho2c, uc, r2beta, delta, r2delta, beta2;
    cout.precision(10);
    std::vector<double> exponent(3), T(10), rho1(10), rho2(10), P(10), u(10), lntc(7), lnrhoc(7), lnuc(7), lnpc(7);
    int stop, w, k, r;
    std::string fline;
    int number_lines = 0;
    w = 0;
    k = 0;
    r = 0;

    cout << "inside beta exponent calculation" << endl;

    //Read env data for beta
    double data[12][5];
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
    convertor >> data[row][col];
    }
    }
    file.close();
/*
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
*/
    //Find critical point to calculate beta
    for(int i=0; i<10; i++)
    {
    //i+1 because first line is header
    T[i] = data[i+1][0];
    rho1[i] = data[i+1][1];
    rho2[i] = data[i+1][2];
    P[i] = data[i+1][3];
    u[i] = data[i+1][4];
    }

    //Tc = data[11][0];
    rho1c = data[11][1];
    rho2c = data[11][2];
    //rhoc = (rho1c+rho2c)/2;
    Pc = data[11][3];
    uc = data[11][4];

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

    (*beta) = linear_regression_angular(lntc, lnrhoc);

    cout << "beta inside = " << beta << endl;
    cout << "r2beta inside = " << r2beta << endl;

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

    cout << "almost out of beta exponents void = " << beta << endl;
}

void beta_exponent2(double *beta, double Tc, double rhoc)
{
    double r2beta, delta, r2delta, beta2;
    cout.precision(10);
    std::vector<double> exponent(3), T(14), rho1(14), rho2(14), P(14), u(14), lntc(14), lnrhoc(14), lnuc(14), lnpc(14);
    int stop, w, k, r;
    std::string fline;
    int number_lines = 0;
    w = 0;
    k = 0;
    r = 0;

    cout << "inside beta exponent calculation" << endl;

    //Read env data for beta
    double data[15][5];
    ifstream file("../Planilhas de análise/env_exponent2.csv");
    cout << "Reading data to calculate critical exponent beta..." << endl;
    for(int row = 0; row < 15; ++row)
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
    convertor >> data[row][col];
    }
    }
    file.close();
/*
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
*/
    //Find critical point to calculate beta
    for(int i=0; i<14; i++)
    {
    //i+1 because first line is header
    T[i] = data[i+1][0];
    rho1[i] = data[i+1][1];
    rho2[i] = data[i+1][2];
    P[i] = data[i+1][3];
    u[i] = data[i+1][4];
    }

    //Tc = data[14][0];
    //rho1c = data[14][1];
    //rho2c = data[14][2];
    //rhoc = (rho1c+rho2c)/2;
    //Pc = data[14][3];
    //uc = data[14][4];

    cout << "beta exponents 2" << endl;
    cout << "Tc = " << Tc << endl;
    cout << "rhoc = " << rhoc << endl;

    //Calculate beta exponent

        //1. Between 0.2% and 1.5% below Critical Temperature

        //2. Between 0.95 and 0.99 Tr

        //3. 1K below Tc to Tc
        for(int i=0; i<14; i++)
        {
        lntc[i] = log(fabs(T[i]-Tc)/Tc);
        lnrhoc[i] = log((rho2[i]-rho1[i])/rhoc);
        }

    (*beta) = linear_regression_angular(lntc, lnrhoc);

    cout << "beta inside 2 = " << beta << endl;
    cout << "r2beta inside 2 = " << r2beta << endl;

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

    cout << "almost out of beta exponents2 void = " << beta << endl;
}

void delta_exponent(double *delta, double rhoc, double uc)
{
    std::vector<double> u(1000), rho(1000), drho(1000);
    std::vector<double> lndrhoc, lnduc;

    //Read env data for beta
    double data[1000][8];
    ifstream file("../Planilhas de análise/crit_isotherm.csv");
    cout << "Reading data to calculate critical exponent delta..." << endl;
    for(int row = 0; row < 1000; ++row)
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
    convertor >> data[row][col];
    }
    }
    file.close();

    for(int i=0; i<1000; i++)
    {
    rho[i] = data[i][0];
    u[i] = data[i][3];
    drho[i] = fabs((rho[i]-rhoc)/rhoc);

        if(drho[i]>=0.1 && drho[i]<=0.3)
        {
        lndrhoc.push_back(log(drho[i]));
        lnduc.push_back(log(fabs((u[i]-uc)/uc)));
        }
    }

    (*delta) = linear_regression_angular(lndrhoc, lnduc);

    cout << "delta inside = " << delta << endl;
}

vector<double> estimate_L_phi(int k, double T)
{
    //variables


    //read experimental data


    //compare actual T with



}

double fugacity_renormalized(int phase, double xint, double Pint, double bm, double R, double T, double **d2p, double **d2u,
                             vector<double> x, vector<double> rho, double **Pmat, double **umat, double *rho_phase, double *Z_phase)
{
    //std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), x(200);
    std::vector<double> Pvec(1000), Pfvec(1000);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double rho_out, uint, Z, ures, lnphi, phi, rhoint;

    //Adjusting parameters to splin2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    //Interpolating isotherm for given x
    double Pvint;

    splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,0.0005,&Pvint);
    Pvec[0] = Pvint;
    Pfvec[0] = Pvec[0] - Pint;

        //cout << "interpolation: " << 0 << " x " << xint << " rho " << 0.0005 << " " << 0.0005/bm
        //     << " P " << Pvec[0] << " Pf " << Pfvec[0] << endl;

    for(int i=1;i<1000;i++)
    {
        rhoint = rho[i];
        rhoint = double(i)/1000;

        splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,rhoint,&Pvint);
        Pvec[i] = Pvint;
        Pfvec[i] = Pvec[i] - Pint;

        //cout << "interpolation: " << i << " x " << xint << " rho " << rhoint << " " << rhoint/bm
        //     << " P " << Pvec[i] << " Pf " << Pfvec[i] << endl;
    }
    //cout << "interpolated phase: " << phase << " x = " << xint << endl;

    //Bracketing the real densities at given P
    int max1 = 0;
    int min1 = 900;
    int max2 = max1+1;
    int min2 = min1-1;
    while(Pfvec[max1]*Pfvec[max2]>0)
    max2 = max2+5;
    if(max2-10<0) max1 = 0;
    else max1 = max2-10;

    while(Pfvec[min1]*Pfvec[min2]>0)
    min2 = min2-5;
    min1 = min2+10;
    //cout << "bracketed: " << min1 << " " << min2 << " " << max1 << " " << max2 << endl;

    //Calculate coexistence densities in interpolated isotherm for given P
    double rho_vap = falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-5);
    double rho_liq = falsi_spline(rho, Pfvec, rho[min1], rho[min2], 1e-5);
    cout << "densities v/l found: " << rho_vap/bm << " " << rho_liq/bm << endl;

    //Select desired density
    if(phase==1) rho_out = rho_liq;
    if(phase==2) rho_out = rho_vap;

    (*rho_phase) = rho_out;

    //Interpolate chemical potential, calculate residual
    //uint = bicubic_s_int(x,rho,xint,rho_out,umat);

    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);

    ures = uint - R*T*log(rho_out/bm); //Use bm because throughout the function rho is adimensional
    //cout << "chemical potential: " << uint << " ures: " << ures << " bm: "  << bm << endl;

    //Calculate Z
    Z = Pint/(R*T*(rho_out/bm));
    (*Z_phase) = Z;
    //cout << "Z: " << Z << " T: " << T << " P: " << Pint << " R: " << R << " rho_out/bm: " << rho_out/bm << endl;

    //Calculate phi
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);

    //cout << phase << "  phi: " << phi << " " << lnphi << endl;
    cout << "V= " << 1/(rho_out/bm) << " Z= " << Z << " phase: " << phase << endl;
    cout << "phi = " << phi << endl;
    //delete[] datap;
    //delete[] datau;
    //delete[] Pmat;
    //delete[] umat;
/*
        for(int i=0; i<201; i++)
        {
        delete[] datap[i];
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] Pmat[i];
        delete[] umat[i];
        }
        delete[] datap;
        delete[] datau;
        delete[] Pmat;
        delete[] umat;
*/
    return phi;
}

void V_renormalized(int phase, double xint, double Pint, double bm, double R, double T, double **d2p, double **d2u,
                    vector<double> x, vector<double> rho, double **Pmat, double **umat, double *V)
{
    //std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), x(200);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double rho_out, uint, Z, ures, lnphi, phi, rhoint;

    //Adjusting parameters to splin2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    std::vector<double> Pvec(n2), Pfvec(n2);

    //Interpolating isotherm for given x
    double Pvint;

    splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,0.0001,&Pvint);
    Pvec[0] = Pvint;
    Pfvec[0] = Pvec[0] - Pint;

        //cout << "interpolation: " << 0 << " x " << xint << " rho " << 0.0005 << " " << 0.0005/bm
        //     << " P " << Pvec[0] << " Pf " << Pfvec[0] << endl;

    for(int i=1;i<n2;i++)
    {
        rhoint = rho[i];
        rhoint = double(i)/n2;

        splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,rhoint,&Pvint);
        Pvec[i] = Pvint;
        Pfvec[i] = Pvec[i] - Pint;

        //cout << "interpolation: " << i << " x " << xint << " rho " << rhoint << " " << rhoint/bm
        //     << " P " << Pvec[i] << " Pf " << Pfvec[i] << endl;
    }
    //cout << "interpolated phase: " << phase << " x = " << xint << endl;

    //Bracketing the real densities at given P
    int max1 = 0;
    int min1 = 0.9*n2;
    int max2 = max1+1;
    int min2 = min1-1;
    while(Pfvec[max1]*Pfvec[max2]>0)
    max2 = max2+int(n2/200);
    if(max2-int(n2/100)<0) max1 = 0;
    else max1 = max2-int(n2/100);

    while(Pfvec[min1]*Pfvec[min2]>0)
    min2 = min2-int(n2/200);
    min1 = min2+int(n2/100);
    //cout << "bracketed: " << min1 << " " << min2 << " " << max1 << " " << max2 << endl;

    //Calculate coexistence densities in interpolated isotherm for given P
    double rho_vap = falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-7);
    double rho_liq = falsi_spline(rho, Pfvec, rho[min1], rho[min2], 1e-7);
    cout << "densities adimensional: " << rho_vap << " " << rho_liq << endl;
    cout << "densities v/l found: " << rho_vap/bm << " " << rho_liq/bm << endl;

    //Select desired density
    if(phase==1) rho_out = rho_liq;
    if(phase==2) rho_out = rho_vap;

    (*V) = 1/(rho_out/bm);
}


double fugacity_renormalized1(int phase, double xint, double Pint, double bm, double R, double T, double **d2p, double **d2u,
                             vector<double> x, vector<double> rho, double **Pmat, double **umat, double rho_out, double Z)
{
    //std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), x(200);
    std::vector<double> Pvec(1000), Pfvec(1000);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double uint, ures, lnphi, phi, rhoint;

    //Adjusting parameters to splin2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }
/*
    //Interpolating isotherm for given x
    double Pvint;

    splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,0.0005,&Pvint);
    Pvec[0] = Pvint;
    Pfvec[0] = Pvec[0] - Pint;

        //cout << "interpolation: " << 0 << " x " << xint << " rho " << 0.0005 << " " << 0.0005/bm
        //     << " P " << Pvec[0] << " Pf " << Pfvec[0] << endl;

    for(int i=1;i<1000;i++)
    {
        rhoint = rho[i];
        rhoint = double(i)/1000;

        splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,rhoint,&Pvint);
        Pvec[i] = Pvint;
        Pfvec[i] = Pvec[i] - Pint;

        //cout << "interpolation: " << i << " x " << xint << " rho " << rhoint << " " << rhoint/bm
        //     << " P " << Pvec[i] << " Pf " << Pfvec[i] << endl;
    }
    //cout << "interpolated phase: " << phase << " x = " << xint << endl;

    //Bracketing the real densities at given P
    int max1 = 0;
    int min1 = 900;
    int max2 = max1+1;
    int min2 = min1-1;
    while(Pfvec[max1]*Pfvec[max2]>0)
    max2 = max2+5;
    if(max2-10<0) max1 = 0;
    else max1 = max2-10;

    while(Pfvec[min1]*Pfvec[min2]>0)
    min2 = min2-5;
    min1 = min2+10;
    //cout << "bracketed: " << min1 << " " << min2 << " " << max1 << " " << max2 << endl;

    //Calculate coexistence densities in interpolated isotherm for given P
    double rho_vap = falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-5);
    double rho_liq = falsi_spline(rho, Pfvec, rho[min1], rho[min2], 1e-5);
    cout << "densities v/l found: " << rho_vap/bm << " " << rho_liq/bm << endl;

    //Select desired density
    if(phase==1) rho_out = rho_liq;
    if(phase==2) rho_out = rho_vap;

    //Interpolate chemical potential, calculate residual
    //uint = bicubic_s_int(x,rho,xint,rho_out,umat);
*/
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);

    ures = uint - R*T*log(rho_out/bm); //Use bm because throughout the function rho is adimensional
    //cout << "chemical potential: " << uint << " ures: " << ures << " bm: "  << bm << endl;

    //Calculate Z
    //Z = Pint/(R*T*(rho_out/bm));
    cout << "Z: " << Z << " T: " << T << " P: " << Pint << " R: " << R << " rho_out/bm: " << rho_out/bm << endl;

    //Calculate phi
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);

    cout << phase << "  phi: " << phi << " " << lnphi << endl;
    cout << "V= " << 1/(rho_out/bm) << " Z= " << Z << " phase: " << phase << endl;
    cout << "phi = " << phi << endl;
    //delete[] datap;
    //delete[] datau;
    //delete[] Pmat;
    //delete[] umat;
/*
        for(int i=0; i<201; i++)
        {
        delete[] datap[i];
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] Pmat[i];
        delete[] umat[i];
        }
        delete[] datap;
        delete[] datau;
        delete[] Pmat;
        delete[] umat;
*/
    return phi;
}

double fugacity_renormalized2(int phase, double xint, double Pint, double bm, double R, double T, double **d2p, double **d2u,
                             vector<double> x, vector<double> rho, double **Pmat, double **umat, double V)
{
    //std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), x(200);
    std::vector<double> Pvec(1000), Pfvec(1000);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double uint, ures, lnphi, phi, rhoint, Z, rho_out;

    //Adjusting parameters to splin2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }
/*
    //Interpolating isotherm for given x
    double Pvint;

    splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,0.0005,&Pvint);
    Pvec[0] = Pvint;
    Pfvec[0] = Pvec[0] - Pint;

        //cout << "interpolation: " << 0 << " x " << xint << " rho " << 0.0005 << " " << 0.0005/bm
        //     << " P " << Pvec[0] << " Pf " << Pfvec[0] << endl;

    for(int i=1;i<1000;i++)
    {
        rhoint = rho[i];
        rhoint = double(i)/1000;

        splin2(x1a,x2a,Pmat,d2p,n1,n2,xint,rhoint,&Pvint);
        Pvec[i] = Pvint;
        Pfvec[i] = Pvec[i] - Pint;

        //cout << "interpolation: " << i << " x " << xint << " rho " << rhoint << " " << rhoint/bm
        //     << " P " << Pvec[i] << " Pf " << Pfvec[i] << endl;
    }
    //cout << "interpolated phase: " << phase << " x = " << xint << endl;

    //Bracketing the real densities at given P
    int max1 = 0;
    int min1 = 900;
    int max2 = max1+1;
    int min2 = min1-1;
    while(Pfvec[max1]*Pfvec[max2]>0)
    max2 = max2+5;
    if(max2-10<0) max1 = 0;
    else max1 = max2-10;

    while(Pfvec[min1]*Pfvec[min2]>0)
    min2 = min2-5;
    min1 = min2+10;
    //cout << "bracketed: " << min1 << " " << min2 << " " << max1 << " " << max2 << endl;

    //Calculate coexistence densities in interpolated isotherm for given P
    double rho_vap = falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-5);
    double rho_liq = falsi_spline(rho, Pfvec, rho[min1], rho[min2], 1e-5);
    cout << "densities v/l found: " << rho_vap/bm << " " << rho_liq/bm << endl;

    //Select desired density
    if(phase==1) rho_out = rho_liq;
    if(phase==2) rho_out = rho_vap;

    //Interpolate chemical potential, calculate residual
    //uint = bicubic_s_int(x,rho,xint,rho_out,umat);
*/
    rho_out = (1/V)*bm;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,xint*rho_out,&uint);

    ures = uint - R*T*log(xint*rho_out/bm); //Use bm because throughout the function rho is adimensional
    cout << "chemical potential: " << uint << " ures: " << ures << " bm: "  << bm << endl;

    //Calculate Z
    Z = Pint/(R*T*(rho_out/bm));
    cout << "Z: " << Z << " T: " << T << " P: " << Pint << " R: " << R << " rho_out/bm: " << rho_out << endl;

    //Calculate phi
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);

    cout << phase << endl;
    cout << "V= " << 1/(rho_out/bm) << " Z= " << Z << " ures: " << ures << " uint: " << uint << " ures/R/T: " << ures/R/T << endl;
    cout << "phi = " << phi << endl;
    //delete[] datap;
    //delete[] datau;
    //delete[] Pmat;
    //delete[] umat;
/*
        for(int i=0; i<201; i++)
        {
        delete[] datap[i];
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] Pmat[i];
        delete[] umat[i];
        }
        delete[] datap;
        delete[] datau;
        delete[] Pmat;
        delete[] umat;
*/
    return phi;
}

double partial_density(double xint, VectorXd b, double a)
{
    std::vector<double> rho(1000), V(1000), x(1000);
    VectorXd xv(1000);
    double xv2[1000];
    double dVdx, Vx, Vint, rhoint;

    cout << "enter thy realm of partial density" << endl;
    cout << "b: " << b(0) << " " << b(1) << " /xint: " << xint << " /a: " << a << endl;

    x[0] = 1e-6;
    for(int i=0; i<1000; i=i+1)
    {
    x[i] = double(i)/1000;
    if(i==0) x[0] = 1e-6;
    rho[i] = x[i]/b(0)+(1-x[i])/b(1);
    V[i] = 1/rho[i];
    cout << "i / x / rho / V: " << i << " " << x[i] << " " << rho[i] << " " << V[i] << endl;
    }

    dVdx = cspline_deriv1(x,V,xint);
    Vx = cspline(x,V,xint);

    Vint = Vx - dVdx*a;
    rhoint = 1/Vint;

    cout << "dVdx / Vx / Vint / rhoint: " << dVdx << " " << Vx << " " << Vint << " " << rhoint << endl;
    cout << "out rhoint-----------------------" << endl;

    return rhoint;
}


double fugacity_renormalized3(int phase, double xint, double Pint, double bm, double R, double T, double **d2p, double **d2u,
                             vector<double> x, vector<double> rho, double **Pmat, double **umat, double V, double b, double ii,
                             VectorXd bvec)
{
    //std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), x(200);
    std::vector<double> Pvec(1000), Pfvec(1000);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double uint, ures, lnphi, phi, rhoint, Z, rho_out;

    //Adjusting parameters to splin2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    //rho_out = (1/V)/b;
    //rho_out = rho_out*xint;
    //rho_out = partial_density(xint,bvec,ii);
/*
    cout << "1. (1/V)*bm \n";
    rho_out = (1/V)*bm;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================" << endl;

    cout << "2. (1/V)*bm*xint \n";
    rho_out = (1/V)*bm*xint;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================" << endl;

    cout << "3. (1/V)*b \n";
    rho_out = (1/V)*b;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================" << endl;

    cout << "4. (1/V)*b*xint \n";
    rho_out = (1/V)*b*xint;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================";

    cout << "5. (1/V)*b1 \n";
    rho_out = (1/V)*bvec(1);
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================" << endl;

    cout << "3. (1/V)*b1*xint \n";
    rho_out = (1/V)*bvec(1)*xint;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);
    ures = uint;
    Z = Pint*V/(R*T);
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);
    cout << "Z: |  ures:  |  rho_out:  |  xint: \n";
    cout << Z << " " << ures << " " << rho_out << " " << xint << endl;
    cout << "phi = \n" << phi << endl;
    cout << "=====================";
    cout << "=====================" << endl;
*/

    //===============================================================
    rho_out = (1/V)*bm;
    //cout << "rho_out: " << rho_out/bm << endl;
    //cout << "rho_out*bm: " << rho_out << endl;
    splin2(x1a,x2a,umat,d2u,n1,n2,xint,rho_out,&uint);

    ures = uint; //Use bm because throughout the function rho is adimensional
    //cout << "chemical potential: " << uint << " ures: " << ures << " bm: "  << bm << endl;

    //Calculate Z
    Z = Pint*V/(R*T);
    //cout << "Z: " << Z << " T: " << T << " P: " << Pint << " R: " << R << " rho_out/bm: " << rho_out/bm << endl;

    //Calculate phi
    lnphi = ures/R/T - log(Z);
    phi = exp(lnphi);

    //cout << "V= " << V << " Z= " << Z << " ures: " << ures << " xint: " << xint << " ures/R/T: " << ures/R/T << " " << phase << endl;
    //cout << "phi = " << phi << endl;
    //cout << "-----\n";
    //delete[] datap;
    //delete[] datau;
    //delete[] Pmat;
    //delete[] umat;
/*
        for(int i=0; i<201; i++)
        {
        delete[] datap[i];
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] Pmat[i];
        delete[] umat[i];
        }
        delete[] datap;
        delete[] datau;
        delete[] Pmat;
        delete[] umat;
*/
    return phi;
}


void d2Pgen(double **d2y, int n)
{
    std::vector<double> rho(n), x(200);
    std::string fline;
    int number_linesp = 0;

    ifstream filep("../Planilhas de análise/xpp.csv");

    while (std::getline(filep, fline))
        ++number_linesp;
        filep.close();

    //Reading Data Bank for p and u
    double **datap;
    datap = new double *[201];
    for(int k = 0; k <201; k++)
        datap[k] = new double[n+1];

    double **Pmat;
    Pmat = new double *[200];
    for(int k = 0; k <200; k++)
        Pmat[k] = new double[n];

    filep.open("../Planilhas de análise/xpp.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(filep, line);
    if ( !filep.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datap[row][col];
    }
    }
    filep.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datap[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datap[0][j+1];
        Pmat[i][j] = datap[i+1][j+1];
        }

    }

    //Adjusting parameters to splie2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    splie2(x1a,x2a,n1,n2,Pmat,d2y);

        for(int i=0; i<201; i++)
        {
        delete[] datap[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] Pmat[i];
        }
        delete[] datap;
        delete[] Pmat;
}

void d2ugen(double **d2y, int n)
{
    std::vector<double> rho(n), x(200);
    std::string fline;
    int number_linesu = 0;

    ifstream fileu("../Planilhas de análise/xpu.csv");

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    double **umat;
    umat = new double *[200];
    for(int k = 0; k <200; k++)
        umat[k] = new double[n];

    fileu.open("../Planilhas de análise/xpu.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datau[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datau[0][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

    //Adjusting parameters to splie2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    splie2(x1a,x2a,n1,n2,umat,d2y);

        for(int i=0; i<201; i++)
        {
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] umat[i];
        }
        delete[] datau;
        delete[] umat;
}

void d2u1gen(double **d2y, int n)
{
    std::vector<double> rho(n), x(200);
    std::string fline;
    int number_linesu = 0;

    ifstream fileu("../Planilhas de análise/xpu1.csv");

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    double **umat;
    umat = new double *[200];
    for(int k = 0; k <200; k++)
        umat[k] = new double[n];

    fileu.open("../Planilhas de análise/xpu1.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datau[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datau[0][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

    //Adjusting parameters to splie2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    splie2(x1a,x2a,n1,n2,umat,d2y);

        for(int i=0; i<201; i++)
        {
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] umat[i];
        }
        delete[] datau;
        delete[] umat;
}

void d2u2gen(double **d2y, int n)
{
    std::vector<double> rho(n), x(200);
    std::string fline;
    int number_linesu = 0;

    ifstream fileu("../Planilhas de análise/xpu2.csv");

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    double **umat;
    umat = new double *[200];
    for(int k = 0; k <200; k++)
        umat[k] = new double[n];

    fileu.open("../Planilhas de análise/xpu2.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datau[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datau[0][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

    //Adjusting parameters to splie2
    int n1 = x.size();
    int n2 = rho.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = rho[i];
    }

    splie2(x1a,x2a,n1,n2,umat,d2y);

        for(int i=0; i<201; i++)
        {
        delete[] datau[i];
        }
        for(int i=0; i<200; i++)
        {
        delete[] umat[i];
        }
        delete[] datau;
        delete[] umat;
}

void renorm_mat_reader(double **Pmat, double **umat, int n)
{
    std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), x(200);
    std::vector<double> Pvec(n), Pfvec(n);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double rho_out, uint, ures, Z, lnphi, phi, rhoint;

    ifstream filep("../Planilhas de análise/xpp.csv");
    ifstream fileu("../Planilhas de análise/xpu.csv");

    while (std::getline(filep, fline))
        ++number_linesp;
        filep.close();

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datap;
    datap = new double *[201];
    for(int k = 0; k <201; k++)
        datap[k] = new double[n+1];

    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    filep.open("../Planilhas de análise/xpp.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(filep, line);
    if ( !filep.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datap[row][col];
    }
    }
    filep.close();

    fileu.open("../Planilhas de análise/xpu.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datap[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datap[0][j+1];
        Pmat[i][j] = datap[i+1][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

}

void renorm_uu_reader(double **u1mat, double **u2mat, int n)
{
    std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), x(200);
    std::vector<double> Pvec(n), Pfvec(n);
    std::string fline;
    int number_linesu1 = 0;
    int number_linesu2 = 0;
    double rho_out, uint, ures, Z, lnphi, phi, rhoint;

    ifstream fileu1("../Planilhas de análise/xpu1.csv");
    ifstream fileu2("../Planilhas de análise/xpu2.csv");

    while (std::getline(fileu1, fline))
        ++number_linesu1;
        fileu1.close();

    while (std::getline(fileu2, fline))
        ++number_linesu2;
        fileu2.close();

    //Reading Data Bank for p and u
    double **datau1;
    datau1 = new double *[201];
    for(int k = 0; k <201; k++)
        datau1[k] = new double[n+1];

    double **datau2;
    datau2 = new double *[201];
    for(int k = 0; k <201; k++)
        datau2[k] = new double[n+1];

    fileu1.open("../Planilhas de análise/xpu1.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu1, line);
    if ( !fileu1.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau1[row][col];
    }
    }
    fileu1.close();

    fileu2.open("../Planilhas de análise/xpu2.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu2, line);
    if ( !fileu2.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau2[row][col];
    }
    }
    fileu2.close();

    for(int i=0; i<200; i++)
    {
        for(int j=0; j<n; j++)
        {
        u1mat[i][j] = datau1[i+1][j+1];
        u2mat[i][j] = datau2[i+1][j+1];
        }
    }

}

vector<double> renorm_xvec(int n)
{
    std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), x(200);
    std::vector<double> Pvec(n), Pfvec(n);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double rho_out, uint, ures, Z, lnphi, phi, rhoint;

    ifstream filep("../Planilhas de análise/xpp.csv");
    ifstream fileu("../Planilhas de análise/xpu.csv");

    while (std::getline(filep, fline))
        ++number_linesp;
        filep.close();

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datap;
    datap = new double *[201];
    for(int k = 0; k <201; k++)
        datap[k] = new double[n+1];

    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    double **Pmat;
    Pmat = new double *[200];
    for(int k = 0; k <200; k++)
        Pmat[k] = new double[n];

    double **umat;
    umat = new double *[200];
    for(int k = 0; k <200; k++)
        umat[k] = new double[n];

    filep.open("../Planilhas de análise/xpp.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(filep, line);
    if ( !filep.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datap[row][col];
    }
    }
    filep.close();

    fileu.open("../Planilhas de análise/xpu.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datap[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datap[0][j+1];
        Pmat[i][j] = datap[i+1][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

    return x;
}

vector<double> renorm_rhovec(int n)
{
    std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), x(200);
    std::vector<double> Pvec(n), Pfvec(n);
    std::string fline;
    int number_linesp = 0;
    int number_linesu = 0;
    double rho_out, uint, ures, Z, lnphi, phi, rhoint;

    ifstream filep("../Planilhas de análise/xpp.csv");
    ifstream fileu("../Planilhas de análise/xpu.csv");

    while (std::getline(filep, fline))
        ++number_linesp;
        filep.close();

    while (std::getline(fileu, fline))
        ++number_linesu;
        fileu.close();

    //Reading Data Bank for p and u
    double **datap;
    datap = new double *[201];
    for(int k = 0; k <201; k++)
        datap[k] = new double[n+1];

    double **datau;
    datau = new double *[201];
    for(int k = 0; k <201; k++)
        datau[k] = new double[n+1];

    double **Pmat;
    Pmat = new double *[200];
    for(int k = 0; k <200; k++)
        Pmat[k] = new double[n];

    double **umat;
    umat = new double *[200];
    for(int k = 0; k <200; k++)
        umat[k] = new double[n];

    filep.open("../Planilhas de análise/xpp.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(filep, line);
    if ( !filep.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datap[row][col];
    }
    }
    filep.close();

    fileu.open("../Planilhas de análise/xpu.csv");
    for(int row = 0; row < 201; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < n+1; ++col)
    {
    string val;
    getline(iss, val, ';');
    if ( !iss )
    break;
    stringstream convertor(val);
    convertor >> datau[row][col];
    }
    }
    fileu.close();

    for(int i=0; i<200; i++)
    {
        x[i] = datap[i+1][0];
        for(int j=0; j<n; j++)
        {
        if(i<1) rho[j] = datap[0][j+1];
        Pmat[i][j] = datap[i+1][j+1];
        umat[i][j] = datau[i+1][j+1];
        }

    }

    return rho;
}

double di_calculation(double **bfnl, double **bfns, double xkni, double dy, int imax, int jmax, int ngrid)
{
      //cout << "dy = " << dy << endl;
      double **argl;
      argl = new double *[ngrid];
        for(int k=0; k<ngrid; k++)
        argl[k] = new double[ngrid];

      double **args;
      args = new double *[ngrid];
        for(int k=0; k<ngrid; k++)
        args[k] = new double[ngrid];

      double aminl, amins, gnl, gns, qnl, qns, al, as, xkn;
      xkn = 1/xkni;

// ****************************************
// *** search for maximum of gnl and ns ***
// ****************************************

      aminl = 0.0;
      amins = 0.0;
      for (int i=0; i<=imax; i++)
      {
         for (int j=0; j<=jmax; j++)
         {
            gnl = 0.5 * (bfnl[imax+i][jmax+j]+bfnl[imax-i][jmax-j])-bfnl[imax][jmax];
            gns = 0.5 * (bfns[imax+i][jmax+j]+bfns[imax-i][jmax-j])-bfns[imax][jmax];
            argl[i][j]=gnl*xkni;
            args[i][j]=gns*xkni;
            if (argl[i][j] < aminl) aminl=argl[i][j];
            if (args[i][j] < amins) amins=args[i][j];
         }
      }
      //cout << "gnl, gns ok" << endl;

// *****************************
// *** numerical integration ***
// *****************************

      qnl=0.25 * (exp(-argl[0][0]+aminl)+exp(-argl[imax][0]+aminl)+exp(-argl[0][jmax]+aminl)+exp(-argl[imax][jmax]+aminl));
      qns=0.25 * (exp(-args[0][0]+amins)+exp(-args[imax][0]+amins)+exp(-args[0][jmax]+amins)+exp(-args[imax][jmax]+amins));

      //cout << "integrating...0" << endl;

      if(imax>1)
      {
      for(int i=1; i<=imax-1; i++)
      {
         al=argl[i][0]-aminl;
         as=args[i][0]-amins;
         if (al < 30.0) qnl=qnl+0.5*exp(-al);
         if (as < 30.0) qns=qns+0.5*exp(-as);
         al=argl[i][jmax]-aminl;
         as=args[i][jmax]-amins;
         if (al < 30.0) qnl=qnl+0.5*exp(-al);
         if (as < 30.0) qns=qns+0.5*exp(-as);
         //cout << "i: " << i << endl;
      }
      }

      //cout << "integrating...1" << endl;

      if(jmax>1)
      {
      for(int j=1; j<=jmax-1; j++)
      {
         al=argl[0][j]-aminl;
         as=args[0][j]-amins;
         if (al < 30.0) qnl=qnl+0.5*exp(-al);
         if (as < 30.0) qns=qns+0.5*exp(-as);
         al=argl[imax][j]-aminl;
         as=args[imax][j]-amins;
         if (al < 30.0) qnl=qnl+0.5*exp(-al);
         if (as < 30.0) qns=qns+0.5*exp(-as);
         //cout << "j: " << j << endl;
      }
      }

      //cout << "integrating...2" << endl;

      if(imax>1 && jmax>1)
      {
      for (int i=1; i<=imax-1; i++)
      {
         for (int j=1; j<=jmax-1; j++)
         {
            al=argl[i][j]-aminl;
            as=args[i][j]-amins;
            if (al < 30.0) qnl=qnl+exp(-al);
            if (as < 30.0) qns=qns+exp(-as);
         }
      }
      }
      //cout << "integrating...3" << endl;

      double di;

      qnl=log(qnl*dy)-aminl;
      qns=log(qns*dy)-amins;
      di=qns-qnl;
      cout << "imax,jmax: " << imax << " " << jmax << endl;
      cout << "aminl,amins,dy: " << aminl << " " << amins << " " << dy << endl;
      cout << "di,qnl,qns: " << di << " " << qnl << " " << qns << endl;
      //cout << "everything ok" << endl;

// ***********************
// *** end executables ***
// ***********************

      return di;
}

void d2_chem_p(double **fn, double *rho1_vec, double *rho2_vec, int n, double drho1, double drho2,
               double **p2, double **fn111, double **fn222)
{
      double tmp;
      int i, j;

      double **p;
      double **fn_mat;
      double **fm_mat;
      double **fn1;
      double **fn11;
      double **fn2;
      double **fn22;

      p = new double *[n];
        for(int k=0; k<n; k++)
        p[k] = new double[n];

      fn_mat = new double *[n];
        for(int k=0; k<n; k++)
        fn_mat[k] = new double[n];

      fm_mat = new double *[n];
        for(int k=0; k<n; k++)
        fm_mat[k] = new double[n];

      fn1 = new double *[n];
        for(int k=0; k<n; k++)
        fn1[k] = new double[n];

      fn11 = new double *[n];
        for(int k=0; k<n; k++)
        fn11[k] = new double[n];

      fn2 = new double *[n];
        for(int k=0; k<n; k++)
        fn2[k] = new double[n];

      fn22 = new double *[n];
        for(int k=0; k<n; k++)
        fn22[k] = new double[n];


      for (i=0; i<n; i++)
      {
         for (j=0; j<n; j++)
         {
            fn_mat[i][j]=fn[i][j];
            fm_mat[j][i]=fn_mat[i][j];
         }
      }

      //'renorm: spline fitting free energy'
      splie2(rho1_vec,rho2_vec,n,n,fn_mat,fn11);
      splie2(rho1_vec,rho2_vec,n,n,fm_mat,fn22);

      for (i=0; i<n-1; i++)
      {
         for (j=i+1; j<n; j++)
         {
            tmp=fn11[i][j];
            fn11[i][j]=fn11[j][i];
            fn11[j][i]=tmp;
         }
      }

// ****************************************************
// *** compute first derivatives of the free energy ***
// ****************************************************

      //write(6,*) 'renorm: computing df/drho1'
      for (j=0; j<n; j++)
      {
         for (i=0; i<n-1; i++)
         {
            fn1[i][j]=(fn[i+1][j]-fn[i][j])/drho1-drho1*(2.0*fn11[i][j]+fn11[i+1][j])/6.0;
         }
         fn1[n][j]=(fn[n][j]-fn[n-1][j])/drho1+drho1*(fn11[n-1][j]+2.0*fn11[n][j])/6.0;
      }

      //write(6,*) 'renorm: computing df/drho2'
      for (i=0; i<n; i++)
      {
         for (j=0; j<n-1; j++)
         {
            fn2[i][j]=(fn[i][j+1]-fn[i][j])/drho2-drho2*(2.0*fn22[i][j]+fn22[i][j+1])/6.0;
         }
         fn2[i][n]=(fn[n][j]-fn[n-1][j])/drho2+drho2*(fn22[i][n-1]+2.0*fn22[i][n])/6.0;
      }

// **************************************
// *** spline fit chemical potentials ***
// **************************************

      splie2(rho1_vec,rho2_vec,n,n,fn1,fn111);
      splie2(rho1_vec,rho2_vec,n,n,fn2,fn222);

// ****************************************
// *** compute pressure and derivatives ***
// ****************************************

      for (i=0; i<n; i++)
      {
         for (j=0; j<n; j++)
         {
            p[i][j]=rho1_vec[i]*fn1[i][j]+rho2_vec[j]*fn2[i][j]-fn[i][j];
         }
      }

// --- spline fit pressure
      splie2(rho1_vec,rho2_vec,n,n,p,p2);

      //Freeing memory
        delete[] p;
        delete[] fn1;
        delete[] fn2;
        delete[] fn11;
        delete[] fn22;

        for(i=0; i<n; i++)
        {
        delete[] p;
        delete[] fn1;
        delete[] fn2;
        delete[] fn11;
        delete[] fn22;
        }

        delete[] p;
        delete[] fn1;
        delete[] fn2;
        delete[] fn11;
        delete[] fn22;

      //end
}

#endif // RENORMALIZATION_H_INCLUDED
