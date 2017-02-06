#ifndef BROYDN_H_INCLUDED
#define BROYDN_H_INCLUDED

#include <iostream>
#include "numerical.h"
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

void envelope_tracer(double tol, double choice)
{
cout.precision(10);
std::vector<double> rho(1000), P(1000), V(1000), A(1000), Pa(1000), f(1000), u(1000), zero(3), zerom(3);
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

    while(w*1000-1<number_lines)
    {
        cout << "BEGIN======================" << endl;

    //Reading Data Bank
    double data[1000][8];
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

    for(int i=0; i<1000; i++)
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
    zero = dens_newt(rho,f,P,u,tol);
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
    minline = w*1000+1;
    maxline = minline+999;
    }
}




























#endif // ENVELOPE_H_INCLUDED
