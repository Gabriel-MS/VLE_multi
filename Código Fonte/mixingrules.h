#ifndef MIXINGRULES_H_INCLUDED
#define MIXINGRULES_H_INCLUDED

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

#include "Gibbs.h"

using namespace Eigen;
using namespace std;

double a_mixing_rules_function(int nc, int MR, VectorXd ai, VectorXd x, double k12, double b, VectorXd bi, double T, VectorXd q, VectorXd r, MatrixXd A, double R, MatrixXd alfa_NRTL, int EdE, int G_ex_model)
{
    int i, j, nc_nc;
    double a, pre_a, C_, G_ex;

    switch(EdE)
    {
        case 1:  C_ = -0.69314; //SRK , ln(2)
        case 2:  C_ = -0.62323; //PR
        case 3:  C_ = -0.69314; //CPA , ln(2)
    }

    nc_nc = pow(nc,nc);
    MatrixXd aij(nc,nc), aiaj(nc,nc), raiz_aiaj(nc,nc), axx(nc,nc), raiz_aiaj_kij(nc,nc), kij(nc,nc), kij1(nc,nc);


    nc_nc = pow(nc,nc);

    VectorXd one(nc_nc);
        for(i=0; i<nc_nc; i++)
    {
         one(i) = 1;
    }

    kij(0,0) = 0;
    kij(1,0) = k12;
    kij(0,1) = k12;
    kij(1,1) = 0;

switch(MR)
{
    case 3: //HV
    switch(G_ex_model)
    {
    case 1: G_ex = gibbs_excess_UNIQUAC(nc, T, q, r, x, A, R); //UNIQUAC
    case 2: G_ex = gibbs_excess_NRTL(nc, T, R, x, A, alfa_NRTL); //NRTL
    }
    pre_a = (ai*(bi.asDiagonal().inverse())).transpose()*x;
    pre_a = x.transpose()*((bi.asDiagonal().inverse())*ai);
    a = b*(pre_a + G_ex/C_);

    //cout << "G_ex = " << G_ex << endl;
    //cout << "pre_a = " << pre_a << endl;
    //cout << "am = " << a << endl;
    break;

    case 1: //VdW1
    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         aij(i,j) = ai(j);
        }
    }
    aiaj = aij.transpose()*ai.transpose();
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    axx = ((raiz_aiaj_kij*x.asDiagonal()).transpose())*(x.asDiagonal());
    VectorXd a_vector(Map<VectorXd>(axx.data(), axx.cols()*axx.rows()));
    a = one.transpose()*a_vector;
    //a = x.transpose()*(raiz_aiaj_kij*x);
    //cout << "kij = " << endl;
    //cout << "kij1 = " << kij1 << endl;
    //cout << "a = " << a << endl;
    break;
    //--------------------------------------------------------------------------------
    /*case 2: //VdW2
    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         aij(i,j) = ai(j);
        }
    }

    aiaj = aij.transpose()*ai.transpose();
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    axx = ((raiz_aiaj_kij*x.asDiagonal()).transpose())*(x.asDiagonal());
    a_vector(Map<VectorXd>(axx.data(), axx.cols()*axx.rows()));
    a = one.transpose()*a_vector;
    break;
    */
    //---------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
}
     return a;
}

double b_mixing_rules_function(int nc, VectorXd bi, VectorXd x, int MR)
{
    int i, j, nc_nc;
    double b;
    VectorXd xb(nc_nc);

    nc_nc = pow(nc,nc);

    VectorXd one(nc), one_2(nc_nc);

    MatrixXd bij(nc,nc), bibj(nc,nc), bxx(nc,nc);


        for(i=0; i<nc; i++)
    {
         one(i) = 1;
    }

        for(i=0; i<nc_nc; i++)
    {
         one_2(i) = 1;
    }


switch(MR)
{
    //-----------------------------------------------------------------
     case 3: //HV
     b = bi.transpose()*x;
     //cout << "bi = " << bi << endl;
     //cout << "bm = " << b << endl;
     break;
    //-----------------------------------------------------------------

    case 1: //VdW1
     //b = bi.transpose()*x;

     for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         bij(i,j) = bi(j);
        }
    }
    bibj = (bij+(bij.transpose())).array()/2;
    bxx = ((bibj*x.asDiagonal()).transpose())*(x.asDiagonal());
    VectorXd b_vector(Map<VectorXd>(bxx.data(), bxx.cols()*bxx.rows()));
    b = one_2.transpose()*b_vector;
     break;

     //---------------------------------------------------------------
/*
     case 2: //VdW2

    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         bij(i,j) = ai(j);
        }
    }
    bibj = bij.transpose()*ai.transpose();
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    axx = ((raiz_aiaj_kij*x.asDiagonal()).transpose())*(x.asDiagonal());
    VectorXd  a_vector(Map<VectorXd>(axx.data(), axx.cols()*axx.rows()));
    a = one.transpose()*a_vector;

    break;*/

}

     return b;
}


#endif // MIXINGRULES_H_INCLUDED
