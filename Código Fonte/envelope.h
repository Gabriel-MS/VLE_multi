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

    ofstream envelope("../Planilhas de análise/Envelope.csv");
    ifstream file("../Planilhas de análise/Renormalization.csv");

    std::vector<double> V(1000), A(1000), Pa(1000), rho(1000), P(1000), f(1000), u(1000), dens(3);
    std::string fline;
    int stop, w, k, searchline;
    int number_lines = 0;
    double T;
    w = 1;
    k = 0;

    envelope << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << endl;

while (std::getline(file, fline)) ++number_lines;
cout << "number of lines in file: " << number_lines << endl;

    //Reading File
/*
    double data[1000][8];
    for(int row = 0; row < number_lines; ++row)
    {
    if(row >= minline || row <= maxline)
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
    }
*/

double minline = 1;
double maxline = 1000;
int r = 0;

    //Tracing envelope
    while(w*1000<number_lines)
    {

    int row = 0;
    double data[1000][8];
    while(row < 61001)
    {
    cout << row << endl;
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

        if(row >= minline && row <= maxline)
        {
        convertor >> data[r][col];
        r++;
        cout << r << endl;
        }
    }
    cout << "not out" << endl;
    ++row;
    }

    r = 0;

    cout << "w: " << w << endl;

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
    cout << "rho = " << rho.size() << " / " << rho[0] << " / " << rho[999] << endl;
    //cout << "P = " << P.size() << " / " << P[999] << endl;
    //cout << "f = " << f.size() << " / " << f[999] << endl;
    //cout << "u = " << u.size() << " / " << u[999] << endl;

std::vector<double> dens(3);
double f1, f2, u1, u2, P1, P2;
if(choice==1) dens = dens_maxwell(rho, P);
if(choice==2)
{
    dens = dens_maxwell(rho, P);
    //dens[0] = dens[0]/bm;
    //dens[1] = dens[1]/bm;
    //dens[2] = dens[2]/bm/bm*am;
}
if(choice==3) dens = dens_area(V, A, P);
if(choice==4) dens = dens_newt(rho, f, P, u, tol);
if(choice==5)
{
    dens = dens_newt(rho, f, P, u, tol);
    //dens[0] = dens[0]/bm;
    //dens[1] = dens[1]/bm;
    //dens[2] = dens[2]/bm/bm*am;
}

    //zero = dens_maxwell(rho, P);
    //cout << "maxwell = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    //zero = dens_area(V, A, P);
    //cout << "area = " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;
    //zero = dens_newt(rho,f,P,u,1e-3);
    //cout << "zero = " << T << " / " << zero[0] << " / " << zero[1] << " / " << zero[2] << endl;

    envelope << T << ";" << dens[0] << ";" << dens[1] << ";" << dens[2] << endl;

    w++;
    minline = w*1000+1;
    maxline = minline+999;
    }

}




























#endif // ENVELOPE_H_INCLUDED
