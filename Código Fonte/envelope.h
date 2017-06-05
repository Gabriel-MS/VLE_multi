#ifndef BROYDN_H_INCLUDED
#define BROYDN_H_INCLUDED

#include <iostream>
#include "numerical.h"
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
/*
void envelope_mix_vle(double tol, int choice)
{
cout.precision(10);
std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), zero(6), zerom(7), x(200);
std::vector<double> zerom_last(3), zero_last(6);
int stop, w, k, r;
std::string fline;
int number_linesp = 0;
int number_linesu = 0;
double T, u1, u2;
w = 0;
k = 0;

ofstream Envelope("../Planilhas de análise/vle.csv");
ifstream filep("../Planilhas de análise/xpp.csv");
ifstream fileu("../Planilhas de análise/xpp.csv");
Envelope << "T" << ";" << "P" << ";" << "x" << ";" << "y" << endl;


while (std::getline(filep, fline))
        ++number_linesp;
        cout << "number of lines in file xpp: " << number_linesp << endl;
        filep.close();

while (std::getline(fileu, fline))
        ++number_linesu;
        cout << "number of lines in file xuu: " << number_linesu << endl;
        fileu.close();


    cout << "BEGIN VLE======================" << endl;

    //Reading Data Bank for p
    double datap[201][1001];
    double datau[201][1001];

    filep.open("../Planilhas de análise/xpp.csv");
    for(int row = 0; row < number_linesp; ++row)
    {
    string line;
    getline(filep, line);
    if ( !filep.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < 1001; ++col)
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
    for(int row = 0; row < number_linesp; ++row)
    {
    string line;
    getline(fileu, line);
    if ( !fileu.good() )
    break;
    stringstream iss(line);
    for (int col = 0; col < 1001; ++col)
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
        for(int j=0; j<n; j++)
        {
        x[i] = datap[i+1][0];
        if(i<1) rhob[i] = datap[0][j+1];
        P[i][j] = datap[i][j];
        u[i][j] = datau[i][j];
        }
    }
    //cout << " rho = " << rho[0] << " / " << rho[500] << " / " << rho[999] << endl;
    //cout << " f = " << f[0] << " / " << f[500] << " / " << f[999] << endl;
    //cout << " P = " << P[0] << " / " << P[500] << " / " << P[999] << endl;
    //cout << " u = " << u[0] << " / " << u[500] << " / " << u[999] << endl;

    switch(choice)
    {
    case 1: //Maxwell
        //zerom = dens_maxwell(rho, P, f, tol);
        //cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
        break;

    case 2: //Area
        break;

    case 3: //Newton
        zero = dens_newt(rho,f,P,u,tol);
        cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 4: //Maxwell + Newton
        //zerom = dens_maxwell(rho, P, f, tol);
        //cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
        //zero = dens_newt(rho,f,P,u,tol);
        //cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 5: //Faster Newton
        /zero = dens_newt5(rho,f,P,u,tol);
        //cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 6: //Newton seed
        /if(w==0) zero = dens_newt(rho,f,P,u,tol);
        //else zero = dens_newt_seed(rho,f,P,u,tol,zero_last);
        //cout << "newton seed= " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;
    }

    for(int j=0; j<4; j++) zero_last[j] = zero[j];

    //zero = dens_area(V, A, P);
    //cout << "area = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    Envelope << T << ";" << zero[0] << ";" << zero[1] << ";" << zero[2] << ";" << zero[4] << ";" << fabs(zero[3]-zero[2])
             << ";" << fabs(zero[4]-zero[5]) << endl;
    Envelope_max << T << ";" << zerom[0] << ";" << zerom[1] << ";" << zerom[2] << ";" << zerom[3] << ";" << zerom[4] << ";"
                 << zerom[5] << ";" << zerom[6] << ";" << fabs(zerom[3]-zerom[4]) << ";" << fabs(zerom[5]-zerom[6]) << ";" << endl;

    w++;
    minline = w*n+1;
    maxline = minline+999;
    }
}
*/

void envelope_tracer_seed(double tol, int choice, int n)
{
cout.precision(10);
std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), zero(3), zerom(3);
std::vector<double> zerom_last(3), zero_last(3);
int stop, w, k, r;
std::string fline;
int number_lines = 0;
double T, u1, u2;
w = 0;
k = 0;

ofstream Envelope("../Planilhas de análise/env.csv");
ofstream Envelope_max("../Planilhas de análise/env_max.csv");
ifstream file("../Planilhas de análise/Renormalization.csv");
Envelope << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << endl;
Envelope_max << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << endl;


while (std::getline(file, fline))
        ++number_lines;
        cout << "number of lines in file: " << number_lines << endl;
int minline = 1;
int maxline = minline+999;
r = 0;
file.close();

    while(w*n-1<number_lines)
    {
        cout << "BEGIN======================" << endl;

    //Reading Data Bank
    double data[n][8];
    file.open("../Planilhas de análise/Renormalization.csv");
    for(int row = 0; row < number_lines; ++row)
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
        if(row>=minline && row<=maxline)
        {
        convertor >> data[r][col];
        }
    }
    if(row>=minline && row<=maxline) r++;
    }
    file.close();
    r = 0;

    for(int i=0; i<n; i++)
    {
    T = data[i][7];
    rho[i] = data[i][0];
    P[i] = data[i][4];
    V[i] = 1/rho[i];
    A[i] = data[i][1]/rho[i];
    f[i] = data[i][1];
    u[i] = data[i][3];
    }

    if(w==0)
    {
    //zerom = dens_maxwell(rho, P, tol);
    //cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
    //zero = dens_area(V, A, P);
    //cout << "area = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    zero = dens_newt(rho,f,P,u,tol,n);
    cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    }

    else
    {
    //zerom = dens_maxwell_seed(rho, P, tol, zerom_last);
    //cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
    //zero = dens_area(V, A, P);
    //cout << "area = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    zero = dens_newt_seed(rho,f,P,u,tol,zero_last);
    cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    }

    for(int j=0; j<4; j++)
    {
        zerom_last[j] = zerom[j];
        zero_last[j] = zero[j];
    }

    Envelope << T << ";" << zero[0] << ";" << zero[1] << ";" << zero[2] << endl;
    Envelope_max << T << ";" << zerom[0] << ";" << zerom[1] << ";" << zerom[2] << endl;

    w++;
    minline = w*n+1;
    maxline = minline+999;
    }
}

void envelope_tracer(double tol, int choice, int n)
{
cout.precision(10);
std::vector<double> rho(n), P(n), V(n), A(n), Pa(n), f(n), u(n), zero(6), zerom(7);
std::vector<double> zerom_last(3), zero_last(6);
int stop, w, k, r;
std::string fline;
int number_lines = 0;
double T, u1, u2;
w = 0;
k = 0;

ofstream Envelope("../Planilhas de análise/env.csv");
ofstream Envelope_max("../Planilhas de análise/env_max.csv");
ifstream file("../Planilhas de análise/Renormalization.csv");
Envelope << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << ";" << "delta_P" << ";" << "delta_u" << endl;
Envelope_max << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "P1" << ";" << "P2" << ";" << "P2" << ";"
             << "u1" << ";" << "u2" << ";" << "delta_P" << ";" << "delta_u" << endl;


while (std::getline(file, fline))
        ++number_lines;
        cout << "number of lines in file: " << number_lines << endl;
int minline = 1;
int maxline = minline+999;
r = 0;
file.close();

    //It was w*n-1
    while(w*n<number_lines)
    {
        cout << "BEGIN======================" << endl;

    //Reading Data Bank
    double data[n][8];
    file.open("../Planilhas de análise/Renormalization.csv");
    for(int row = 0; row < number_lines; ++row)
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
        if(row>=minline && row<=maxline)
        {
        convertor >> data[r][col];
        }
    }
    if(row>=minline && row<=maxline) r++;
    }
    file.close();
    r = 0;

    for(int i=0; i<n; i++)
    {
    T = data[i][7];
    rho[i] = data[i][0];
    P[i] = data[i][4];
    V[i] = 1/rho[i];
    A[i] = data[i][1]/rho[i];
    f[i] = data[i][1];
    u[i] = data[i][3];
    }
    //cout << " rho = " << rho[0] << " / " << rho[500] << " / " << rho[999] << endl;
    //cout << " f = " << f[0] << " / " << f[500] << " / " << f[999] << endl;
    //cout << " P = " << P[0] << " / " << P[500] << " / " << P[999] << endl;
    //cout << " u = " << u[0] << " / " << u[500] << " / " << u[999] << endl;

    switch(choice)
    {
    case 1: //Maxwell
        zerom = dens_maxwell(rho, P, f, tol);
        cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
        break;

    case 2: //Area
        break;

    case 3: //Newton
        zero = dens_newt(rho,f,P,u,tol,n);
        cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 4: //Maxwell + Newton
        zerom = dens_maxwell(rho, P, f, tol);
        cout << "maxwell = " << T << " / " << zerom[0] << " / " << zerom[1] << " / " << zerom[2] << endl;
        zero = dens_newt(rho,f,P,u,tol,n);
        cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 5: //Faster Newton
        zero = dens_newt5(rho,f,P,u,tol);
        cout << "newton = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;

    case 6: //Newton seed
        if(w==0) zero = dens_newt(rho,f,P,u,tol,n);
        else zero = dens_newt_seed(rho,f,P,u,tol,zero_last);
        cout << "newton seed= " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
        break;
    }

    for(int j=0; j<4; j++) zero_last[j] = zero[j];

    //zero = dens_area(V, A, P);
    //cout << "area = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    Envelope << T << ";" << zero[0] << ";" << zero[1] << ";" << zero[2] << ";" << zero[4] << ";" << fabs(zero[3]-zero[2])
             << ";" << fabs(zero[4]-zero[5]) << endl;
    Envelope_max << T << ";" << zerom[0] << ";" << zerom[1] << ";" << zerom[2] << ";" << zerom[3] << ";" << zerom[4] << ";"
                 << zerom[5] << ";" << zerom[6] << ";" << fabs(zerom[3]-zerom[4]) << ";" << fabs(zerom[5]-zerom[6]) << ";" << endl;

    w++;
    minline = w*n+1;
    maxline = minline+999;
    }
}

vector<double> T_tracer(double T, double dnew, int flag, double rho)
{
    std::vector<double> out(2);

    if(flag==2)
    {
    T = T + 0.01;

    if(isnan(rho)==1) flag=5;
    if(isnan(rho)==1) T = T - 0.02;
    }

    if(flag==1)
    {

        if(isnan(rho)==1)
        {
            T = T - 0.09;
            flag = 2;
        }

        else
        {
            T = T + 0.1;
        }

    }

    if(flag==0)
    {
    if(isnan(dnew)==1)
    {
        if(isnan(dnew)==1)
        {
            T = T - 1.9;
            flag = 1;
        }

        else
        {
            T = T + 0.1;
        }
    }

    else
    {
        T = T + 2;
    }

    }

    out[0] = T;
    out[1] = flag;
    return out;
}

vector<double> T_tracer_CPA(double T, double dnew, int flag, double rho)
{
    std::vector<double> out(2);

    if(flag==3)
    {
    T = T + 0.001;

    if(isnan(rho)==1) flag=5;
    if(isnan(rho)==1) T = T - 0.002;
    }

    if(flag==2)
    {
        if(isnan(rho)==1)
        {
            T = T - 0.009;
            flag = 3;
        }

        else
        {
            T = T + 0.01;
        }
    }

    if(flag==1)
    {

        if(isnan(rho)==1)
        {
            T = T - 0.09;
            flag = 2;
        }

        else
        {
            T = T + 0.1;
        }

    }

    if(flag==0)
    {
    if(isnan(dnew)==1)
    {
        if(isnan(dnew)==1)
        {
            T = T - 1.9;
            flag = 1;
        }

        else
        {
            T = T + 0.1;
        }
    }

    else
    {
        T = T + 2;
    }

    }

    out[0] = T;
    out[1] = flag;
    return out;
}

vector<double> T_tracer_inflexion(int EdE, double T, int flag)
{
    std::vector<double> out(2);
    double Temp;

    switch(flag)
    {
        case 0:
        Temp = T - 0.9;
        flag = 1;
            break;

        case 1:
        Temp = T - 0.09;
        flag = 2;
            break;

        case 2:
        if(EdE==3)
        {
            Temp = T - 0.009;
            flag = 3;
        }

        else
        {
            Temp = T - 0.01;
            flag = 6;
        }
            break;

        case 3:
        Temp = T - 0.001;
        flag = 6;
            break;

    }

    out[0] = Temp;
    out[1] = flag;

    return out;
}


#endif // ENVELOPE_H_INCLUDED
