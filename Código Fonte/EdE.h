#ifndef EDE_H_INCLUDED
#define EDE_H_INCLUDED


#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

#include "Gibbs.h"
#include "mixingrules.h"
#include "Association.h"

using namespace Eigen;
using namespace std;


//Function to calculate alfa
VectorXd alfa_function(int EdE, int nc, VectorXd Tr, VectorXd omega, VectorXd a0, VectorXd c1)
{
    VectorXd alfa(nc), pre_Tr(nc), alfa_termo1(nc), omega2(nc), pre_alfa2(nc);
    MatrixXd pre_alfa(nc,nc);

switch(EdE)
	{
    case 1: //SRK
        pre_Tr = -Tr.array().pow(0.5)+1;
        omega2 = omega.array().pow(2);
        alfa_termo1 = 0.48 + omega.array()*1.574 - 0.176*omega2.array();
        pre_alfa = (alfa_termo1 * pre_Tr.transpose());
        pre_alfa2 = 1 + pre_alfa.diagonal().array();
        //pre_alfa = 1 + ((omega.array()*1.574 - 0.176*omega.array().pow(2) + 0.480) * pre_Tr.asDiagonal()).diagonal();
        alfa = pre_alfa2.array().pow(2);
        break;

    case 2: //PR
        pre_Tr = 1-Tr.array().pow(0.5);
        omega2 = omega.array().pow(2);
        alfa_termo1 = 0.37464 + omega.array()*1.54226 - 0.26992*omega2.array();
        pre_alfa = (alfa_termo1 * pre_Tr.transpose());
        pre_alfa2 = 1 + pre_alfa.diagonal().array();
        //pre_alfa = 1 + ((0.37464 + omega*1.54226 - omega.array().pow(2)*0.26992) * pre_Tr.asDiagonal()).diagonal();
        alfa = pre_alfa2.array().pow(2);
        break;

    case 3: //CPA
        pre_Tr = (-Tr.array().pow(0.5))+1;
        omega2 = omega.array().pow(2);
        alfa_termo1 = c1;
        pre_alfa = (alfa_termo1 * pre_Tr.transpose());
        pre_alfa2 = 1 + pre_alfa.diagonal().array();
        //pre_alfa = 1 + ((omega.array()*1.574 - 0.176*omega.array().pow(2) + 0.480) * pre_Tr.asDiagonal()).diagonal();
        alfa = pre_alfa2.array().pow(2);
        break;
	}
return alfa;
}

//Function to get set up EdE parameters
VectorXd EdE_parameters_function(int EdE)
{
    double sigma, epsilon, OMEGA, PSI;
    VectorXd EdE_parameters(4);

    switch(EdE)
	{
    case 1: //SRK
        sigma = 1;
        epsilon = 0;
        OMEGA = 0.08664;
        PSI = 0.42748;
        break;

    case 2: //PR
        sigma = 1 + pow(2, 0.5);
        epsilon = 1 - pow(2, 0.5);
        OMEGA = 0.07780;
        PSI = 0.45724;
        break;

    case 3: //CPA-SRK
        sigma = 1;
        epsilon = 0;
        OMEGA = 0.08664;
        PSI = 0.42748;
        break;
	}
EdE_parameters[0] = sigma;
EdE_parameters[1] = epsilon;
EdE_parameters[2] = OMEGA;
EdE_parameters[3] = PSI;

	return EdE_parameters;
}

//Function to calculate ai
VectorXd a_function(int nc, double R, VectorXd EdE_parameters, VectorXd omega, VectorXd Tc, VectorXd Pc, VectorXd alfa, int EdE, VectorXd a0)
{
        double sigma, epsilon, OMEGA, PSI;
        VectorXd pre_a(nc), a(nc), Tc2(nc), pre_a2(nc);
        //MatrixXd pre_a2(nc,nc);

switch(EdE)
{
    case 1: //SRK
    sigma = EdE_parameters[0];
    epsilon = EdE_parameters[1];
    OMEGA = EdE_parameters[2];
    PSI = EdE_parameters[3];

    pre_a = PSI*R*R*alfa.array();
    Tc2 = Tc.array().pow(2)/1000; //sem essa divisão, o valor extrapola o limite superior
    pre_a2 = pre_a.transpose()*Tc2.asDiagonal();
    a = pre_a2.transpose()*Pc.asDiagonal().inverse();
    a = a.array()*1000;

    pre_a = PSI*R*R*alfa.array();
    Tc2 = Tc.array().pow(2)/1000; //sem essa divisão, o valor extrapola o limite superior
    pre_a2 = (Tc2.asDiagonal())*pre_a;
    a = (Pc.asDiagonal().inverse())*pre_a2;
    a = a.array()*1000;
    break;

    case 2: //PR
    sigma = EdE_parameters[0];
    epsilon = EdE_parameters[1];
    OMEGA = EdE_parameters[2];
    PSI = EdE_parameters[3];

    pre_a = PSI*R*R*alfa.array();
    Tc2 = Tc.array().pow(2)/1000; //sem essa divisão, o valor extrapola o limite superior
    pre_a2 = pre_a.transpose()*Tc2.asDiagonal();
    a = pre_a2.transpose()*Pc.asDiagonal().inverse();

    pre_a = PSI*R*R*alfa.array();
    Tc2 = Tc.array().pow(2)/1000; //sem essa divisão, o valor extrapola o limite superior
    pre_a2 = (Tc2.asDiagonal())*pre_a;
    a = (Pc.asDiagonal().inverse())*pre_a2;

    a = a.array()*1000;
    break;

    case 3: //CPA
    a = (a0.asDiagonal())*alfa;
    break;
}

    return a;
}

//Function to calculate bi
VectorXd b_function(int nc, double R, VectorXd EdE_parameters, VectorXd omega, VectorXd Tc, VectorXd Pc, VectorXd alfa, VectorXd bCPA, int EdE)
{
        double sigma, epsilon, OMEGA, PSI;
        VectorXd b(nc), pre_b(nc);

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

switch(EdE)
{
case 1: //SRK
    pre_b = OMEGA*R*Tc.array();
    b = pre_b.transpose()*Pc.asDiagonal().inverse();
    break;

case 2: //PR
    pre_b = OMEGA*R*Tc.array();
    b = pre_b.transpose()*Pc.asDiagonal().inverse();
    break;

case 3: //CPA
    b = bCPA;
    break;
}

    return b;
}

//Objective function to calculate V using CPA-SRK EoS using Michelsen's method (2006)
double CPA_volume_obj_function(int nc, double V, double R, double P, double am, double bm,
                               double T, VectorXd x, VectorXd X, double B, double iota, int iter, VectorXd a)
{
    //Coloquei VectorXd a nas variáveis daqui e de vol function!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!****************************************

    VectorXd one_4(4), one_nc(nc), one_4nc(4*nc);
    MatrixXd pre_F(nc,4);
    double F, pre_F_v;
    int i;

    int j;
    double dlng_dv, P_CPA, eta, g, fV, D_T;
    MatrixXd kij(nc,nc), kij1(nc,nc), aij(nc,nc), aiaj(nc,nc), raiz_aiaj(nc,nc), raiz_aiaj_kij(nc,nc);

    V = bm/iota;

    one_4 <<    1,
                1,
                1,
                1;

    for(i=0; i<nc; i++)
    {
         one_nc(i) = 1;
    }

    for(i=0; i<(4*nc); i++)
    {
         one_4nc(i) = 1;
    }

    pre_F = (x*one_4.transpose()).transpose();
    VectorXd pre_F_vector(Map<VectorXd>(pre_F.data(), pre_F.cols()*pre_F.rows()));
    pre_F_v = pre_F_vector.transpose()*(one_4nc-X);
    F = ((R*T/(V-bm) - am/(V*(V+bm)) - 0.5*R*T/V * (1+0.475*B/(V-0.475*B)) * pre_F_v) - P);
    /*
    cout << "(R*T/(V-bm) - am/(V*(V+bm)) = " << R*T/(V-bm) - am/(V*(V+bm)) << endl;
    cout << "0.5*R*T/V =                   " << 0.5*R*T/V << endl;
    cout << "(1+0.475*B/(V-0.475*B)) =     " << (1+0.475*B/(V-0.475*B)) << endl;
    cout << "pre_F_v =                     " << pre_F_v << endl;
    */
/*
    //falta kij!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //***************************************************************************************************************

    kij << 0, 0,
           0, 0;

    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         aij(i,j) = a(j);
        }
    }
    aiaj = aij.transpose()*a.transpose();
    aiaj = a*(a.transpose());
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    aij = raiz_aiaj_kij;

    dlng_dv = 1/g * (-0.475*B*((1/(V-0.475*B))*(1/(V-0.475*B))));
    eta = bm/(4*V);
    g = 1/(1-1.9*eta);
    fV = -(1/(R*V*(V+B)));
    D_T = x.transpose()*(aij*x);


    P_CPA = R*T/V - R*T*(-B/(V*(V-B)-D_T/T*fV)) - R*T*0.5/V*(1-V*dlng_dv)*pre_F_v;
    F = P - P_CPA;
*/
    F = (1-iota)*(F);
    return F;
}

//Function to calculate volume using CPA-SRK EoS using Michelsen's method (2006)
VectorXd volume_function(int nc, int EdE, int phase, VectorXd x, VectorXd X, VectorXd EdE_parameters, double bm,
                       double am, double R, double T, double P, double tolV, double tolZ, VectorXd b, int combining_rule,
                       MatrixXd beta_row, MatrixXd beta_col, MatrixXd E_row, MatrixXd E_col, VectorXd alfa, double tolX,
                       VectorXd n_v, double *V, double Vinit, VectorXd a, double *V_obj, double *Q, double *dP_dV_output,
                       double BETCR)
{
    int i;
    double B, sigma, epsilon, OMEGA, PSI, Vi, errorV, errorZ, condV, Vnew, F_obj, F_obj_plus, F_obj_minus;
    double Zi, Z, qe;
    int d;
    double V1, k, alpha, beta, gama, f_v, f_v_der;
    MatrixXd pre_F(nc,4);

    tolZ = tolV;

B = bm*P/(R*T);
qe = (am/(bm*R*T));

    VectorXd one_4(4);

    one_4 <<    1,
                1,
                1,
                1;

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    switch(EdE)
    {
    case 1: //SRK

        d = 0;
    switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------
    Vi = R*T/P; //Initial guess
    errorV = tolV + 1;
      while(errorV>tolV)
      {

       alpha = (sigma+epsilon-1)*bm-R*T/P;
       beta = sigma*epsilon*bm*bm-(R*T/P+bm)*(sigma+epsilon)*bm+am/P;
       gama = -(R*T/P+bm)*sigma*epsilon*bm*bm-am*bm/P;

       f_v = 2*Vi*Vi*Vi+alpha*Vi*Vi-gama;
       f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       //V = Vi - f_v/f_v_der;
       (*V) = f_v/f_v_der;
       errorV = fabs((*V)-Vi);
       Vi = (*V);
       if(d==1000)
       {
           errorV = tolV-1;
       }
       d++;
      }
      V1 = (*V);
      k = -(alpha+(*V))/2;
      (*V) = k - pow(k*k+2*k*(*V)-beta,0.5);
      (*V) = min(V1,(*V));
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Vi = R*T/P;
    errorV = tolV + 1;
      while(errorV>tolV)
      {

       alpha = (sigma+epsilon-1)*bm-R*T/P;
       beta = sigma*epsilon*bm*bm-(R*T/P+bm)*(sigma+epsilon)*bm+am/P;
       gama = -(R*T/P+bm)*sigma*epsilon*bm*bm-am*bm/P;

       f_v = 2*Vi*Vi*Vi+alpha*Vi*Vi-gama;
       f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       (*V) = f_v/f_v_der;
       errorV = fabs((*V)-Vi);
       Vi = (*V);
       if(d==1000)
       {
           errorV = tolV-1;
       }
       d++;
      }
      V1 = (*V);
      k = -(alpha+(*V))/2;
      (*V) = k + pow(k*k+2*k*(*V)-beta,0.5);
      (*V) = max(V1,(*V));
    break;
    }
        break;

    case 2: //PR
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Zi = B; //Chute inicial pra fase líquida
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        (*V) = R*T*Z/P;
        break;


        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Zi = 1; //Chute inicial pra fase vapor
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        (*V) = R*T*Z/P;
        break;
        }

    case 3: //CPA
    double iota, iota2;
    double F_obj, F_obj_minus, F_obj_plus, F_obj_derivate, iota_min, iota_max, cond_iota, F_obj_old, Bcpa;
    double iota_old, iota_new, deltaV;
    int d, k;
    long double div;
    VectorXd one_4nc(4*nc);

    double QVV, h, v_dlng_dv, dlng_dv, v_d2lng_dv2, dP_dV_SRK, max_dP_dV;
    VectorXd X2(4*nc), pre_Xnew(4*nc), QXV(4*nc), dX_dV(4*nc), dg_dV(4*nc); //, dPassoc_dV(4*nc), dP_dV(4*nc);

    MatrixXd h_1(4*nc,4*nc), DELTA(4*nc,4*nc), xXD(4*nc,4*nc), K(4*nc,4*nc), QXX1(4*nc,4*nc), QXX(4*nc,4*nc), one_4_x(4*nc,4*nc);

    double dPassoc_dV, dP_dV;

    double dh_dV, d2lng_dv2;
    double Q_func;

    double V_residue, Residue, iota_res;
    VectorXd X_residue(4*nc);

    for(i=0;i<(4*nc);i++)
    {
        one_4nc(i) = 1;
    }
    //*******************************************
    n_v = x;

    k = 0;

    Bcpa = n_v.transpose()*b;
    //Chute inicial de V parte de SRK
    //chutes iniciais para iota

    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        iota = 0.99;
        iota = bm/Vinit;
        //cout << "\n LIQUID" << endl;
        //cout << "Vinit = " << Vinit << endl;
        break;

        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        iota = bm/(bm+(R*T/P));
        iota = bm/Vinit;
        //cout << "\n VAPOR" << endl;
        //cout << "Vinit = " << Vinit << endl;
        break;
        }

        (*V) = bm/iota;

    //Cálculo de V a partir do chute inicial
    iota_min = 0;
    iota_max = 1;

    i=0;
    cond_iota = tolV+1;
    //X = fraction_nbs_initial(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
    //               EdE_parameters, b, tolZ, (*V));
    deltaV = 0;

/*
        ofstream Residue_liquid("../Planilhas de análise/F_obj_liquid.csv");
        ofstream Residue_vapor("../Planilhas de análise/F_obj_vapor.csv");

for (iota_res=0.0001 ; iota_res<=1.000 ; iota_res=iota_res+0.0002)
{
    V_residue = bm/iota_res;

    X_residue = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                     EdE_parameters, b, tolZ, V_residue, deltaV, X, i, a, &Q_func, BETCR);

    Residue = CPA_volume_obj_function(nc, V_residue, R, P, am, bm, T, x, X_residue, Bcpa, iota, i, a);

if(phase==1)
{
Residue_liquid << iota_res << ";" << Residue << endl;
}

if(phase==2)
{
Residue_vapor << iota_res << ";" << Residue << endl;
}
}
//cout << "calculated \n";
//cin.get();
*/


int iter;
iter = 0;
    //while(cond_iota>tolV || max_dP_dV>0)
    //while(cond_iota>tolV || max_dP_dV>0 || i<3)
    while(cond_iota>tolV)
    {

    if(i==0)
    {
    X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                     EdE_parameters, b, tolZ, (*V), deltaV, X, i, a, &Q_func, BETCR);
    }

    F_obj = CPA_volume_obj_function(nc, (*V), R, P, am, bm, T, x, X, Bcpa, iota, i, a);
    F_obj_plus = CPA_volume_obj_function(nc, (*V), R, P, am, bm, T, x, X, Bcpa, iota+iota*0.00001, i, a);
    F_obj_minus = CPA_volume_obj_function(nc, (*V), R, P, am, bm, T, x, X, Bcpa, iota-iota*0.00001, i, a);
    F_obj_derivate =(F_obj_plus-F_obj_minus)/(2*iota*0.00001);
    //F_obj_derivate = (F_obj-F_obj_old)/(iota-iota_old);


    if(i==0)
    {
        F_obj_derivate = 100000000000;
        //F_obj_derivate = F_obj/iota;
    }

    //F_obj_derivate = ((F_obj+0.001)-(F_obj-0.001))/0.002;


    iota_old = iota;
    F_obj_old = F_obj;

    if(F_obj>0)
    {
        iota_max = iota;
    }

    else
    {
        iota_min = iota;
    }
    //div = (i<3?0.001:1.)*F_obj/F_obj_derivate;

    div = F_obj/F_obj_derivate;
/*
    if(i<3)
    {
        div = div*0.00001;
    }
*/
    iota_new = iota - div;

    if(iota_min<iota_new && iota_new<iota_max)
    {
        iota2 = iota_new;
        //cout << "caso 1 \n" << endl;
    }

    else
    {
        iota2 = (iota_min+iota_max)/2;
        //cout << "caso 2 \n" << endl;
    }

    cond_iota = fabs(iota2 - iota)/iota;
    cond_iota = fabs(F_obj);
    //esse if como comentário

    if(i>0)
    {
    cond_iota = fabs(F_obj/F_obj_derivate);
    }

    i = i+1;

    iota = iota2;

    deltaV = bm/iota - (*V);

    (*V) = bm/iota;
    (*V_obj) = F_obj;

    //-----------------------------------------------------------------------------------
    one_4_x = one_4*x.transpose();
    VectorXd one_4_x_vector(Map<VectorXd>(one_4_x.data(), one_4_x.cols()*one_4_x.rows()));

    if(i!=0)
    {
    X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                     EdE_parameters, b, tolZ, (*V), deltaV, X, i, a, &Q_func, BETCR);
    }

    v_dlng_dv = 0.475*bm/((*V)-0.475*bm);
    dlng_dv = v_dlng_dv/(*V);
    v_d2lng_dv2 = 0.475*bm*(2-0.475*bm/(*V))/((*V)*(*V)-0.475*(*V));
    d2lng_dv2 = -0.475*bm/((*V)*((*V)-0.475*bm));
    //********************

    dlng_dv = -0.475*bm/((*V)*((*V)-0.475*bm));
    v_dlng_dv = (*V)*dlng_dv;
    v_d2lng_dv2 = 0.475*bm*(2*(*V)-0.475*bm)/(((*V)*(*V)-0.475*bm*(*V))*((*V)*(*V)-0.475*bm*(*V)));
    d2lng_dv2 = v_d2lng_dv2/(*V);

//********************************ESSA PARTE É NOVA PARA CALCULAR NOVO h!!!!!********************************************************************
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, b, tolZ, (*V), BETCR);
    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/(*V);
    QXV = (one_4_x_vector.asDiagonal()*pre_Xnew)*(1/(*V)-1);
    X2 = X.array().pow(2);
    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/(*V);
    QXX1 = (((X2.asDiagonal().inverse())*one_4_x_vector).asDiagonal());
    QXX = -QXX1-K;
    dX_dV = (QXX.inverse())*(-QXV);
//**********************************************************************************************************************

    h_1 = (n_v*one_4.transpose()).transpose();
    VectorXd h_1_vector(Map<VectorXd>(h_1.data(), h_1.cols()*h_1.rows()));
    h = (h_1_vector.transpose())*(one_4nc-X);

    dh_dV = (h_1_vector.transpose())*(one_4nc-dX_dV);

    QVV = (-h/(2*(*V)))*((1-v_dlng_dv)-(dlng_dv + v_d2lng_dv2));
    //*******************
    QVV = (-h/(2*(*V)*(*V)))*((1-v_dlng_dv)-(dlng_dv + v_d2lng_dv2));
    //QVV utilizando dh_dV
    QVV = (dh_dV*(1/(*V)-dlng_dv) + h*(-1/((*V)*(*V))-d2lng_dv2))/2;

//***************************ANTERIORMENTE ESTA PARTE ESTAVA ATIVA***************************************************************
/*
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, b, tolZ, (*V));
    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/(*V);
    QXV = (one_4_x_vector.asDiagonal()*pre_Xnew)*(1/(*V)-1);
    X2 = X.array().pow(2);
    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/(*V);
    QXX1 = (((X2.asDiagonal().inverse())*one_4_x_vector).asDiagonal());
    QXX = -QXX1-K;
    dX_dV = (-QXX.inverse())*QXV;
*/
//*************************************************************************************************************************

    dg_dV = QXV;
    //dPassoc_dV = (-1/(R*T))*(X.array()*QVV+dg_dV.transpose()*dX_dV);
    //**********************************************************************


    dPassoc_dV = (-R*T)*(QVV+dg_dV.transpose()*dX_dV);


    //**********************************************************************
    dP_dV_SRK = -R*T/(((*V)-bm)*((*V)-bm))+am/(bm*(*V)*(*V))-am/(bm*((*V)+bm)*((*V)+bm));
    dP_dV_SRK = -R*T/(((*V)-bm)*((*V)-bm))+am/(((*V)+bm)*((*V)+bm)*(*V))+am/((*V)*((*V)+bm)*((*V)+bm));
    //dP_dV = dPassoc_dV.array() + dP_dV_SRK;

    //max_dP_dV = dP_dV.maxCoeff();
    max_dP_dV = dP_dV_SRK - (R*T)*QVV;

    //ESSES ESTAVAM ATIVOS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dP_dV = dPassoc_dV + dP_dV_SRK;
    max_dP_dV = dP_dV;

    //MUDEI AQUI################################################################################################3
    /*
    F_obj_plus = CPA_volume_obj_function(nc, (*V), R, P, am, bm, T, x, X, Bcpa, iota+iota*0.00001, i, a);
    F_obj_minus = CPA_volume_obj_function(nc, (*V), R, P, am, bm, T, x, X, Bcpa, iota-iota*0.00001, i, a);
    F_obj_derivate =(F_obj_plus-F_obj_minus)/(2*iota*0.00001);

    max_dP_dV = ((((F_obj_plus)/(1-iota+iota*0.00001))+P)-(((F_obj_minus)/(1-iota-iota*0.00001))+P))/(2*(*V)*0.00001);
    //max_dP_dV = F_obj_derivate;
*/

    (*dP_dV_output) = max_dP_dV;


    if(isnan(*V)==1 || isinf(*V)==1)
    {
        (*V) = Vinit+Vinit*0.001*k;
        iota = bm/(*V);
        X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                     EdE_parameters, b, tolZ, (*V), deltaV, X, i, a, &Q_func, BETCR);
        max_dP_dV = 0.1;
        cond_iota = tolV+1;
        k++;
        cout << "VOLUME = NAN or INF \n";

        Vinit = (*V);
        //cin.get();
    }

    if (i>500)
    {
        max_dP_dV = -0.1;
        cond_iota = tolV-1;
        (*V) = Vinit;
        cout << "VOLUME MAX ITER REACHED \n";
        //cin.get();
    }
    //cout << "F_obj = " << F_obj << endl;
    //cout << "X = \n" << X << endl;
    //cout << "max_dP_dV = " << max_dP_dV << endl;
    //cout << "cond_iota = " << cond_iota << endl;
    //cout << "i = " << i << endl;
    //cin.get();
    }
    //cout << "iota = " << iota << endl;
    //cout << "dP_dV = " << max_dP_dV << endl;
    //cin >> i;
    break;



    }

return X;
}

/*
//Function to calculate volume using CPA-SRK EoS using Michelsen's method (2006)
VectorXd total_volume_function(int nc, int EdE, int phase, VectorXd x, VectorXd X, VectorXd EdE_parameters, double bm,
                       double am, double R, double T, double P, double tolV, double tolZ, VectorXd b, int combining_rule,
                       MatrixXd beta_row, MatrixXd beta_col, MatrixXd E_row, MatrixXd E_col, VectorXd alfa, double tolX,
                       VectorXd n_v, double *V, double Vinit)
{
    int i, d;
    double B, sigma, epsilon, OMEGA, PSI, Vi, errorV, errorZ, condV, Vnew, F_obj, F_obj_plus, F_obj_minus;
    double Zi, Z, qe, iota, iota2, iota_old, iota_new, deltaV;
    double V1, k, alpha, beta, gama, f_v, f_v_der, QVV, h, v_dlng_dv, dlng_dv, v_d2lng_dv2, dP_dV_SRK, max_dP_dV;
    double F_obj, F_obj_minus, F_obj_plus, F_obj_derivate, iota_min, iota_max, cond_iota, F_obj_old, Bcpa;
    MatrixXd pre_F(nc,4), h_1(4*nc,4*nc), DELTA(4*nc,4*nc), xXD(4*nc,4*nc), K(4*nc,4*nc), QXX1(4*nc,4*nc), QXX(4*nc,4*nc), one_4_x(4*nc,4*nc);
    VectorXd one_4nc(4*nc), X2(4*nc), pre_Xnew(4*nc), QXV(4*nc), dX_dV(4*nc), dg_dV(4*nc), dPassoc_dV(4*nc), dP_dV(4*nc);

    tolZ = tolV;

B = bm*P/(R*T);
qe = (am/(bm*R*T));

    VectorXd one_4(4);

    one_4 <<    1,
                1,
                1,
                1;

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    for(i=0;i<(4*nc);i++)
    {
        one_4nc(i) = 1;
    }

    Bcpa = n_v.transpose()*b;
    //Chute inicial de V parte de SRK
    //chutes iniciais para iota
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        iota = 0.99;
        iota = bm/Vinit;
        break;

        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        iota = bm/(bm+(R*T/P));
        iota = bm/Vinit;
        break;
        }

        (*V) = bm/iota;

    //Cálculo de V a partir do chute inicial
    iota_min = 0;
    iota_max = 1;

    i=0;
    deltaV = 0;
    cond_iota = tolV+1;
    X = fraction_nbs_initial(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                    EdE_parameters, b, tolZ, (*V));
int iter;
iter = 0;
    while(cond_iota>tolV || max_dP_dV>0)
    {

    if(i!=0)
    {
    X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row, tolX, x, EdE,
                     EdE_parameters, b, tolZ, (*V), deltaV, X);
    }

    F_obj = ;
    F_obj_derivate = (F_obj-F_obj_old)/(iota-iota_old);

    if(i==0)
    {
        F_obj_derivate = 100000000000;
    }

    F_obj_old = F_obj;

    //-----------------------------------------------------------------------------------
    one_4_x = one_4*x.transpose();
    VectorXd one_4_x_vector(Map<VectorXd>(one_4_x.data(), one_4_x.cols()*one_4_x.rows()));
    v_dlng_dv = 0.475*bm/((*V)-0.475*bm);
    dlng_dv = v_dlng_dv/(*V);
    v_d2lng_dv2 = 0.475*bm*(2-0.475*bm/(*V))/((*V)*(*V)-0.475*(*V));

    h_1 = (n_v*one_4.transpose()).transpose();
    VectorXd h_1_vector(Map<VectorXd>(h_1.data(), h_1.cols()*h_1.rows()));
    h = (h_1_vector.transpose())*(one_4nc-X);

    QVV = (-h/(2*(*V)))*((1-v_dlng_dv)-(dlng_dv + v_d2lng_dv2));

    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, b, tolZ, (*V));
    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/(*V);
    QXV = (one_4_x_vector.asDiagonal()*pre_Xnew)*(1/(*V)-1);
    X2 = X.array().pow(2);
    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/(*V);
    QXX1 = (((X2.asDiagonal().inverse())*one_4_x_vector).asDiagonal());
    QXX = -QXX1-K;
    dX_dV = (-QXX.inverse())*QXV;
    dg_dV = QXV;
    dPassoc_dV = (-1/(R*T))*(X.array()*QVV+dg_dV.transpose()*dX_dV);
    dP_dV_SRK = -R*T/(((*V)-bm)*((*V)-bm))+am/(bm*(*V)*(*V))-am/(bm*((*V)+bm)*((*V)+bm));
    dP_dV = dPassoc_dV.array() + dP_dV_SRK;

    max_dP_dV = dP_dV.maxCoeff();
    max_dP_dV = dP_dV_SRK - (R*T)*QVV;

    if (i>100)
    {
        max_dP_dV = -0.1;
        cond_iota = tolV-1;
        (*V) = Vinit;
    }
    }
    break;


return X;
}

*/

//Function to calculate fugacity
VectorXd fugacity_function(int nc, int phase, double am, double bm, VectorXd a, VectorXd b, double R, double T, double P,
                           double tolZ, VectorXd EdE_parameters, int MR, VectorXd q_prime, VectorXd r, MatrixXd A, VectorXd x,
                           VectorXd qUNIQUAC, int EdE, MatrixXd alfa_NRTL, int G_ex_model, double k12, VectorXd X, double tolV,
                           double V, VectorXd n_v, double Vt, double *Z_phase, double *u_phase, VectorXd y)
{
    //Variables-----------------------------------------------------------------------
    int d;
    double B, q, sigma, epsilon, OMEGA, PSI, Zi, Z, errorZ, errorV, I, logZB, Z1, C_, qe, k, V1;
    double a1, a2, a3, Q, R1, Q3R2, A1, theta, Za, Zb, Zc, Vi, f_v, f_v_der, alpha, beta, gama;
    VectorXd q_(nc), phi(nc), q_I(nc), pre_phi(nc), bZ1(nc), gamma(nc), ln_gamma(nc);
    MatrixXd pre_q_(nc,nc), diag_a(nc,nc);
    //--------------------------------------------------------------------------------

    //C* to use Huron-Vidal mixing rule
    switch(EdE)
    {
        case 1:  C_ = -0.69314; //SRK
        case 2:  C_ = -0.62323; //PR
    }

    //Main part of calculation
    switch(EdE)
    {
    case 1: //SRK
    B = bm*P/(R*T);
    qe = (am/(bm*R*T));

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

d = 0;

    switch(MR) //Calculating mixing rule
    {
    case 1: //VdW
    diag_a = a.asDiagonal();
    pre_q_ = (diag_a/am).array().pow(0.5);
    q_ = qe*(b/bm-2*pre_q_.diagonal());
    break;

/*
    case 2: //VdW2
    break;
*/
    case 3: //HV
    switch(G_ex_model)
    {
        case 1: gamma = gamma_UNIQUAC(nc, qUNIQUAC, q_prime, r, x, A); //UNIQUAC
        case 2: gamma = gamma_NRTL(nc, T, R, x, A, alfa_NRTL);
    }
    ln_gamma = gamma.array().log();
    q_ = -((b.asDiagonal().inverse()*a)/(R*T) + ln_gamma/C_);
    break;
    }

/*
    A1 = am*P/((R*T)*(R*T));
    a1 = -1+B;
    a2 = A1-3*B*B-2*B;
    a3 = -A1*B+B*B+B*B*B;
    Q = (a1*a1-3*a2)/9;
    R1 = (2*a1*a1*a1-9*a1*a2+27*a3)/54;
    Q3R2 = Q*Q*Q-R1*R1;
    theta = acos(R1/(pow(Q,1.5)));

    Za = -2*(pow(Q,0.5))*cos(theta/3)-a1/3;
    Zb = -2*(pow(Q,0.5))*cos((theta+2*3.14159)/3)-a1/3;
    Zc = -2*(pow(Q,0.5))*cos((theta+4*3.14159)/3)-a1/3;

    if(phase == 1)
    {
        Z = min(min(Za,Zb),Zc);
    }

    if(phase == 2)
    {
        Z = max(max(Za,Zb),Zc);
    }

    if(isnan(Z)==1)
    {
       switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------
    Zi = B; //Initial guess
    errorZ = tolZ + 1;
      while(errorZ>tolZ)
      {
       Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
       errorZ = fabs(Z-Zi);
       Zi = Z;
       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Zi = 1; //Initial guess
    errorZ = tolZ + 1;
     while(errorZ>tolZ)
      {
       Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
       errorZ = fabs(Z-Zi);
       Zi = Z;

       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;
    }
    }*/

    switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------

    Vi = R*T/P; //Initial guess
    errorV = tolV + 1;
      while(errorV>tolV)
      {
       //f_v = P-(R*T/(Vi-bm) - am/(Vi*(Vi+bm)));
       //f_v_der = (R*T/((Vi-bm)*(Vi-bm)) - am*(2*Vi+bm)/((Vi*(Vi+bm))*(Vi*(Vi+bm))));

       alpha = (sigma+epsilon-1)*bm-R*T/P;
       beta = sigma*epsilon*bm*bm-(R*T/P+bm)*(sigma+epsilon)*bm+am/P;
       gama = -(R*T/P+bm)*sigma*epsilon*bm*bm-am*bm/P;

       //f_v = Vi*Vi*Vi+alpha*Vi*Vi+beta*Vi+gama;
       //f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       f_v = 2*Vi*Vi*Vi+alpha*Vi*Vi-gama;
       f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       //V = Vi - f_v/f_v_der;
       V = f_v/f_v_der;
       errorV = fabs(V-Vi);
       Vi = V;
       if(d==1000)
       {
           errorV = tolV-1;
       }
       d++;
      }
      V1 = V;
      //V = -(V+alpha)-pow((V+alpha)*(V+alpha)-4*(V*V+alpha*V+beta),0.5);
      //V = V/2;
      k = -(alpha+V)/2;
      V = k - pow(k*k+2*k*V-beta,0.5);
      V = min(V1,V);
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Vi = R*T/P;
    errorV = tolV + 1;
      while(errorV>tolV)
      {
       //f_v = P-(R*T/(Vi-bm) - am/(Vi*(Vi+bm)));
       //f_v_der = (R*T/((Vi-bm)*(Vi-bm)) - am*(2*Vi+bm)/((Vi*(Vi+bm))*(Vi*(Vi+bm))));

       alpha = (sigma+epsilon-1)*bm-R*T/P;
       beta = sigma*epsilon*bm*bm-(R*T/P+bm)*(sigma+epsilon)*bm+am/P;
       gama = -(R*T/P+bm)*sigma*epsilon*bm*bm-am*bm/P;

       //f_v = Vi*Vi*Vi+alpha*Vi*Vi+beta*Vi+gama;
       //f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       f_v = 2*Vi*Vi*Vi+alpha*Vi*Vi-gama;
       f_v_der = 3*Vi*Vi+2*alpha*Vi+beta;

       //V = Vi - f_v/f_v_der;
       V = f_v/f_v_der;
       errorV = fabs(V-Vi);
       Vi = V;
       if(d==1000)
       {
           errorV = tolV-1;
       }
       d++;
      }
      V1 = V;
      //V = -(V+alpha)+pow((V+alpha)*(V+alpha)-4*(V*V+alpha*V+beta),0.5);
      //V = V/2;
      k = -(alpha+V)/2;
      V = k + pow(k*k+2*k*V-beta,0.5);
      V = max(V1,V);
    break;
    }

    Z = P*V/(R*T);
/*
    switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------
    Zi = B; //Initial guess
    errorZ = tolZ + 1;
      while(errorZ>tolZ)
      {
       Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
       errorZ = fabs(Z-Zi);
       Zi = Z;
       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Zi = 1; //Initial guess
    errorZ = tolZ + 1;
     while(errorZ>tolZ)
      {
       Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
       errorZ = fabs(Z-Zi);
       Zi = Z;

       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;
    }

*/
    (*Z_phase) = Z;

    I = (1/(sigma-epsilon))*(log((Z+B*sigma)/(Z+B*epsilon)));
    logZB = log(Z-B);
    q_I = q_*I;
    Z1 = Z-1;
    bZ1 = b/bm*Z1;
    pre_phi = bZ1.array()-logZB+q_I.array();
    phi =  (pre_phi).array().exp(); //Aqui deve-se diagonalizar o vetor phi para uma matriz "A" e depois fazer e^A
    break;

    //==============================================================================================================
    case 2: //PR
    B = bm*P/(R*T);
    qe = (am/(bm*R*T));

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

d = 0;
/*
    switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------
    Zi = B; //Initial guess
    errorZ = tolZ + 1;
      while(errorZ>tolZ)
      {
       Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
       errorZ = fabs(Z-Zi);
       Zi = Z;
       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Zi = 1; //Initial guess
    errorZ = tolZ + 1;
     while(errorZ>tolZ)
      {
       Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
       errorZ = fabs(Z-Zi);
       Zi = Z;

       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;
    }
*/
    switch(MR) //Calculating mixing rule
    {
    case 1: //VdW
    diag_a = a.asDiagonal();
    pre_q_ = (diag_a/am).array().pow(0.5);
    q_ = qe*(b/bm-2*pre_q_.diagonal());
    break;

/*
    case 2: //VdW2
    break;
*/
    case 3: //HV
    switch(G_ex_model)
    {
        case 1: gamma = gamma_UNIQUAC(nc, qUNIQUAC, q_prime, r, x, A); //UNIQUAC
        case 2: gamma = gamma_NRTL(nc, T, R, x, A, alfa_NRTL);
    }
    ln_gamma = gamma.array().log();
    q_ = -((b.asDiagonal().inverse()*a)/(R*T) + ln_gamma/C_);
    break;

    }


    A1 = am*P/((R*T)*(R*T));
    a1 = -1+B;
    a2 = A1-3*B*B-2*B;
    a3 = -A1*B+B*B+B*B*B;
    Q = (a1*a1-3*a2)/9;
    R1 = (2*a1*a1*a1-9*a1*a2+27*a3)/54;
    Q3R2 = Q*Q*Q-R1*R1;
    theta = acos(R1/(pow(Q,1.5)));

    Za = -2*(pow(Q,0.5))*cos(theta/3)-a1/3;
    Zb = -2*(pow(Q,0.5))*cos((theta+2*3.14159)/3)-a1/3;
    Zc = -2*(pow(Q,0.5))*cos((theta+4*3.14159)/3)-a1/3;
    if(phase == 1)
    {
        Z = min(min(Za,Zb),Zc);
    }

    if(phase == 2)
    {
        Z = max(max(Za,Zb),Zc);
    }
    /*
    cout << "fraction = " << x << endl;
    cout << "P = " << P << endl;
    cout << "R = " << R << endl;
    cout << "am = " << am << endl;
    cout << "A1 = " << A1 << endl;
    cout << "a1 = " << a1 << endl;
    cout << "a2 = " << a2 << endl;
    cout << "a3 = " << a3 << endl;
    cout << "Q = " << Q << endl;
    cout << "R1 = " << R1 << endl;
    cout << "Q3R2 = " << Q3R2 << endl;
    cout << "theta = " << theta << endl;

    cout << "Z1 = " << Za << endl;
    cout << "Z2 = " << Zb << endl;
    cout << "Z3 = " << Zc << endl;
*/

    //if(isnan(Z)==1)
    //{
       switch(phase) //Calculating compressibility factor
    {
    case 1: //Liquid phase - PR - SRK -------
    Zi = B; //Initial guess
    errorZ = tolZ + 1;
      while(errorZ>tolZ)
      {
       Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
       errorZ = fabs(Z-Zi);
       Zi = Z;
       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;


    case 2: //Vapor phase - PR - SRK ----------------------------
    Zi = 1; //Initial guess
    errorZ = tolZ + 1;
     while(errorZ>tolZ)
      {
       Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
       errorZ = fabs(Z-Zi);
       Zi = Z;

       if(d==1000)
       {
           errorZ = tolZ-1;
       }
       d++;
      }
    break;
    }
    //}

    (*Z_phase) = Z;

    I = (1/(sigma-epsilon))*(log((Z+B*sigma)/(Z+B*epsilon)));
    logZB = log(Z-B);
    q_I = q_*I;
    Z1 = Z-1;
    bZ1 = b/bm*Z1;
    pre_phi = bZ1.array()-logZB+q_I.array();
    phi =  (pre_phi).array().exp(); //Aqui deve-se diagonalizar o vetor phi para uma matriz "A" e depois fazer e^A
    //cout << "Z = " << Z << endl;
    //cout << "phi = " << phi << endl;
    break;
    //==============================================================================================================

    case 3: //CPA
    int i, j, nc4;
    nc4 = 4*nc;
    double Z_assoc, Z_SRK, Z_CPA, p_dlng_dp, h, Fn, f, fV, D_T, FB, FD, Bcpa, eta, g, u_assoc_4, n_t;
    VectorXd u_assoc(nc), u_CPA(nc), u_SRK(nc), ln_phi(nc), u_assoc_1(nc), lnX(nc), u_assoc_3b(nc), u_assoc_3c(nc), u_assoc_3(4*nc);
    VectorXd one_4(4), one_nc(nc), one_4nc(nc4), Biv(nc), Biv1(nc), Div(nc), dlng_dn(nc);
    MatrixXd one_4c(nc4,nc), kij(nc,nc), aij(nc,nc), aiaj(nc,nc), raiz_aiaj(nc,nc), raiz_aiaj_kij(nc,nc), kij1(nc,nc), bij(nc,nc), h_1(nc,4), u_assoc_2(nc,4);

    double Am, Bm, ln_phi_phys2;
    VectorXd ln_phi_phys1(nc), ln_phi_phys3(nc);


sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    B = bm*P/(R*T);
    qe = (am/(bm*R*T));

    one_4 << 1,
             1,
             1,
             1;

    for(i=0; i<nc; i++)
    {
         one_nc(i) = 1;
    }

    for(i=0; i<(4*nc); i++)
    {
         one_4nc(i) = 1;
    }

    one_4c << 1, 0,
              1, 0,
              1, 0,
              1, 0,
              0, 1,
              0, 1,
              0, 1,
              0, 1;

    if(nc==3)
    {
        one_4c << 1, 0, 0,
                  1, 0, 0,
                  1, 0, 0,
                  1, 0, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 0, 1,
                  0, 0, 1,
                  0, 0, 1,
                  0, 0, 1;
    }

    kij(0,0) = 0;
    kij(1,0) = k12;
    kij(0,1) = k12;
    kij(1,1) = 0;

    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         aij(i,j) = a(j);
        }
    }
    aiaj = aij.transpose()*a.transpose();
    aiaj = a*(a.transpose());
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    aij = raiz_aiaj_kij;
    /*
    cout << "a = \n" << a << endl;
    cout << "aij = \n" << aij << endl;
    cout << "aiaj = \n" << aiaj << endl;
    cout << "raiz_aiaj = \n" << raiz_aiaj << endl;
    */

    //*********************************************************************************************************
    //n_v = x;
    //cout << "n_v phi" << n_v << endl;


    //n_v = n_v.asDiagonal()*x;
    n_v = x;
    double n_t_phase;
    VectorXd z(nc), n_phase(nc);

    z = (x.array()+y.array())/2;
    n_t = one_nc.transpose()*n_v;
    n_phase = (((y.array()-z.array()).abs())/(((y.array()-x.array())).abs()))*n_t;
    n_t_phase = one_nc.transpose()*n_phase;

    //cout << "n_v = " << n_v << endl;
    //cout << "WOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOW" << endl;

    //cout << "n_v after = " << n_v << endl;
    //*********************************************************************************************************
    //cout << "x = \n" << x << endl;
    //cout << "n_v = \n" << n_v << endl;

     //Parâmetros de auxílio para os cálculos de Z e u
     bij = (((b*(one_nc.transpose()))+((b*(one_nc.transpose())).transpose())).array())/2;
     Bcpa = n_v.transpose()*b;


     Bcpa = x.transpose()*b;
     //Bcpa = n_phase.transpose()*b;


     Biv1 = 2*((bij*n_v).array());
     //Biv1 = 2*((bij*n_phase).array());//-----------------NEW-------------------NEW----------------NEW------------------------

     Biv = (Biv1.array()-Bcpa)/n_t;

     Div = 2*((aij*n_v).array());
     //Div = 2*((aij*n_phase).array());//-----------------NEW-------------------NEW----------------NEW------------------------

     Vt = V;

     dlng_dn = (b.array())*(0.475/(V-0.475*Bcpa));


     //dlng_dn = (Biv.array())*(0.475/(V-0.475*Bcpa));


     p_dlng_dp = 0.475*bm/(V-0.475*bm);
     //p_dlng_dp = 0.475*Bcpa/(V-0.475*Bcpa);
     eta = bm/(4*V);
     g = 1/(1-1.9*eta);
     h_1 = (n_v*one_4.transpose()).transpose();

     h_1 = (x*one_4.transpose()).transpose();

     VectorXd h_1_vector(Map<VectorXd>(h_1.data(), h_1.cols()*h_1.rows()));
     h = (h_1_vector.transpose())*(one_4nc-X);

     //Associative compressibility factor
     Z_assoc = -0.5*(1+p_dlng_dp)*h;

     //Physical(SRK) compressibility factor
     Z_SRK = V/(V-bm)-am/(R*T*(V+bm));

     //CPA compressibility factor
     Z_CPA = Z_SRK + Z_assoc;
     (*Z_phase) = Z_CPA;

     //Associative residual chemical potential for CPA
     //g = (1-0.5*eta)/((1-eta)*(1-eta)*(1-eta));
     //dlng_dn = b*P/(R*T)*(2.5-eta)/(8*g*Z_CPA)/((1-eta)*(1-eta)*(1-eta)*(1-eta));


     lnX = X.array().log();
     u_assoc_1 = one_4nc.transpose()*(((lnX.asDiagonal())*one_4c));
     u_assoc_2 = ((dlng_dn)*one_4.transpose()).transpose();
     VectorXd u_assoc_2_vector(Map<VectorXd>(u_assoc_2.data(), u_assoc_2.cols()*u_assoc_2.rows()));
     //u_assoc_3 = u_assoc_2_vector.transpose()*(one_4nc-X); ORIGINAL!!!!!!
     u_assoc_3 = u_assoc_2_vector.asDiagonal()*(one_4nc-X);
     u_assoc_3b = one_4nc.transpose()*(((u_assoc_3.asDiagonal())*one_4c));
     u_assoc_3c = n_v.asDiagonal()*u_assoc_3b;

     //u_assoc_3c = n_phase.asDiagonal()*u_assoc_3b;//-----------------NEW-------------------NEW----------------NEW------------------------


     u_assoc_4 = (one_nc.transpose())*u_assoc_3c;
     u_assoc = u_assoc_1.array()-0.5*u_assoc_4;


     lnX = X.array().log();
     u_assoc_1 = one_4c.transpose()*lnX;
     u_assoc_2 = h*dlng_dn.array();
     u_assoc = u_assoc_1-0.5*u_assoc_2;

     //u_assoc = u_assoc_1.array()+u_assoc_4;
/*
     cout << "x = \n" << x << endl;
     cout << "X = \n" << X << endl;
     cout << "lnX = \n" << lnX << endl;
     cout << "dlng_dn = \n" << dlng_dn << endl;
     cout << "0.5*((one_nc.transpose()*u_assoc_4).array()) = \n" << 0.5*((one_nc.transpose()*u_assoc_4).array()) << endl;
     cout << "u_assoc_3 = \n" << u_assoc_3 << endl;
     cout << "u_assoc_2 = \n" << u_assoc_2 << endl;
     cout << "u_assoc_1.array()-0.5*((one_nc.transpose()*u_assoc_4).array()) = \n" << u_assoc_1.array()-0.5*((one_nc.transpose()*u_assoc_4).array()) << endl;
     cout << "u_assoc_1.array() = \n" << u_assoc_1.array() << endl;
     cout << "u_assoc_2_vector = \n" << u_assoc_2_vector << endl;
     cout << "u_assoc = \n" << u_assoc << endl;
     cin.get();
*/
     //Physical(SRK) residual chemical potential for CPA
     Fn = -(log(1-Bcpa/Vt));
     f = (log(1+Bcpa/Vt))/(R*Bcpa);
     //f = 1/(R*Bcpa)*(log(1+Bcpa/Vt));
     fV = -(1/(R*Vt*(Vt+Bcpa)));
     D_T = n_phase.transpose()*(aij*n_phase);

     D_T = am;

     FB = n_t/(Vt-Bcpa)+D_T*(f+Vt*fV)/(T*Bcpa);
     FD = -f/T;
     u_SRK = Fn + FB*(Biv.array()) + FD*(Div.array());

     //Cálculo do potencial químico residual da CPA
/*
    qe = (am/(bm*R*T));
    diag_a = a.asDiagonal();
    pre_q_ = (diag_a/am).array().pow(0.5);
    q_ = qe*(b/bm-2*pre_q_.diagonal());

    B = bm*P/(R*T);

    I = (1/(sigma-epsilon))*(log((Z_SRK+B*sigma)/(Z_SRK+B*epsilon)));
    logZB = log(Z_SRK-B);
    q_I = q_*I;
    Z1 = Z_SRK-1;
    bZ1 = b/bm*Z1;
    pre_phi = bZ1.array()-logZB+q_I.array();
    u_SRK = pre_phi.array() + log(Z_SRK);
*/
     u_CPA = u_assoc + u_SRK;
     (*u_phase) = u_CPA(0);

     /*
     cout << "V = " << V << endl;

     cout << "Bcpa = " << Bcpa << endl;
     cout << "Biv = " << Biv << endl;
     cout << "Div = " << Div << endl;
     cout << "FB = " << FB << endl;
     cout << "FD = " << FD << endl;
     cout << "D_T = " << D_T << endl;
     cout << "fV = " << fV << endl;
     cout << "f = " << f << endl;
     cout << "Fn = " << Fn << endl;

     cout << "Z_SRK = " << Z_SRK << endl;
     cout << "Z_assoc = " << Z_assoc << endl;
     cout << "Z_CPA = " << Z_CPA << endl;
     cout << "log Z_CPA = " << log(Z_CPA) << endl;
     cout << "u_assoc = " << u_assoc << endl;
     cout << "u_SRK = " << u_SRK << endl;
     */

     ln_phi = u_CPA.array() - log(Z_CPA);
     phi = ln_phi.array().exp();

     Am = am*P/(R*R*T*T);
     Bm = bm*P/(R*T);
     ln_phi_phys1 = b/bm*(Z_CPA-1);
     ln_phi_phys2 = log(Z_CPA-Bm);
     ln_phi_phys3 = 2*(x.transpose()*aij).transpose().array()/am - (b/bm).array();
     ln_phi_phys3 = Am/Bm*ln_phi_phys3*log((Z_CPA+Bm)/Z_CPA);
     ln_phi = ln_phi_phys1.array()-ln_phi_phys2-ln_phi_phys3.array();
     ln_phi = ln_phi.array()+u_assoc.array();
     //phi = ln_phi.array().exp();
     /*
     cout << "///////////////////////////////" << endl;
     cout << "Z_CPA = " << Z_CPA << endl;
     cout << "Am = " << Am << endl;
     cout << "Bm = " << Bm << endl;
     cout << "u_assoc = " << u_assoc << endl;
     cout << "ln_phi_phys1 = " << ln_phi_phys1.transpose() << endl;
     cout << "ln_phi_phys2 = " << ln_phi_phys2 << endl;
     cout << "aij = \n" << aij << endl;
     cout << "am = " << am << endl;
     cout << "x = \n" << n_v << endl;
     cout << "x.transpose()*aij = \n" << x.transpose()*aij << endl;
     cout << "2*(x.transpose()*aij) = \n" << 2*(x.transpose()*aij) << endl;
     cout << "(x.transpose()*aij).transpose().array()/am = \n" << (x.transpose()*aij).transpose().array()/am << endl;
     cout << "2*(x.transpose()*aij).transpose().array()/am = \n" << 2*(x.transpose()*aij).transpose().array()/am << endl;
     cout << "b/bm = \n" << b/bm << endl;
     cout << "ln_phi_phys3 = " << ln_phi_phys3.transpose() << endl;
     cout << "ln_phi = " << ln_phi.transpose() << endl;
     cout << "phi = " << phi.transpose() << endl;
     */
     break;

    }


return phi;
}

//Function to calculate fugacity
VectorXd u_CPA_function(int nc, int phase, double am, double bm, VectorXd a, VectorXd b, double R, double T, double P, double tolZ,
                           VectorXd EdE_parameters, int MR, VectorXd q_prime, VectorXd r, MatrixXd A, VectorXd x, VectorXd qUNIQUAC,
                           int EdE, MatrixXd alfa_NRTL, int G_ex_model, double k12, VectorXd X, double tolV, double V, VectorXd n_v, double Vt)
{
    //Variables-----------------------------------------------------------------------
    double B, q, sigma, epsilon, OMEGA, PSI, Zi, Z, errorZ, I, logZB, Z1, C_, qe;
    VectorXd q_(nc), phi(nc), q_I(nc), pre_phi(nc), bZ1(nc), gamma(nc), ln_gamma(nc);
    MatrixXd pre_q_(nc,nc), diag_a(nc,nc);
    //--------------------------------------------------------------------------------

    //Main part of calculation
    int i, j, nc4;
    nc4 = 4*nc;
    double Z_assoc, Z_SRK, Z_CPA, n_t, p_dlng_dp, h, Fn, f, fV, D_T, FB, FD, Bcpa, eta, g, u_assoc_3, u_assoc_4;
    VectorXd u_assoc(nc), u_CPA(nc), u_SRK(nc), ln_phi(nc), u_assoc_1(nc), lnX(nc), u_assoc_3b(nc);
    VectorXd one_4(4), one_nc(nc), one_4nc(nc4), Biv(nc), Biv1(nc), Div(nc), dlng_dn(nc);
    MatrixXd one_4c(nc4,nc), kij(nc,nc), aij(nc,nc), aiaj(nc,nc), raiz_aiaj(nc,nc), raiz_aiaj_kij(nc,nc), kij1(nc,nc), bij(nc,nc), h_1(nc,4), u_assoc_2(nc,4);

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    B = bm*P/(R*T);
    qe = (am/(bm*R*T));

    one_4 << 1,
             1,
             1,
             1;

    for(i=0; i<nc; i++)
    {
         one_nc(i) = 1;
    }

    for(i=0; i<(4*nc); i++)
    {
         one_4nc(i) = 1;
    }

    one_4c << 1, 0,
              1, 0,
              1, 0,
              1, 0,
              0, 1,
              0, 1,
              0, 1,
              0, 1;

    if(nc==3)
    {
        one_4c << 1, 0, 0,
                  1, 0, 0,
                  1, 0, 0,
                  1, 0, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 1, 0,
                  0, 0, 1,
                  0, 0, 1,
                  0, 0, 1,
                  0, 0, 1;
    }

    kij(0,0) = 0;
    kij(1,0) = k12;
    kij(0,1) = k12;
    kij(1,1) = 0;

    for(i=0; i<nc; i++)
    {
        for(j=0; j<nc; j++)
        {
         aij(i,j) = a(j);
        }
    }
    aiaj = aij.transpose()*a.transpose();
    aiaj = a*(a.transpose());
    raiz_aiaj = aiaj.array().pow(0.5);
    kij1 = -kij.array()+1;
    raiz_aiaj_kij = raiz_aiaj.cwiseProduct(kij1);
    aij = raiz_aiaj_kij;
    /*
    cout << "a = \n" << a << endl;
    cout << "aij = \n" << aij << endl;
    cout << "aiaj = \n" << aiaj << endl;
    cout << "raiz_aiaj = \n" << raiz_aiaj << endl;
    */

    //*********************************************************************************************************
    //n_v = x;
    //cout << "n_v phi" << n_v << endl;
    n_v = n_v.asDiagonal()*x;
    //cout << "n_v after = " << n_v << endl;
    //*********************************************************************************************************
    //cout << "x = \n" << x << endl;
    //cout << "n_v = \n" << n_v << endl;

     //Parâmetros de auílio para os cálculos de Z e u
     n_t = one_nc.transpose()*n_v;
     bij = (((b*(one_nc.transpose()))+((b*(one_nc.transpose())).transpose())).array())/2;
     Bcpa = n_v.transpose()*b;
     Biv1 = 2*((bij*n_v).array());
     Biv = (Biv1.array()-Bcpa)/n_t;
     Div = 2*((aij*n_v).array());
     dlng_dn = (b.array())*(0.475/(Vt-0.475*Bcpa));
     p_dlng_dp = 0.475*bm/(V-0.475*bm);
     p_dlng_dp = 0.475*Bcpa/(V-0.475*Bcpa);
     eta = bm/(4*V);
     g = 1/(1-1.9*eta);
     h_1 = (n_v*one_4.transpose()).transpose();
     VectorXd h_1_vector(Map<VectorXd>(h_1.data(), h_1.cols()*h_1.rows()));
     h = (h_1_vector.transpose())*(one_4nc-X);
     //cout << "h = " << h << endl;

     //Associative compressibility factor
     Z_assoc = -0.5*(1+p_dlng_dp)*h;
     //Z_assoc = +2*(1+eta/g*p_dlng_dp)*h;

     //Physical(SRK) compressibility factor
     Z_SRK = V/(V-bm)-am/(R*T*(V+bm));

     //CPA compressibility factor
     Z_CPA = Z_SRK + Z_assoc;

     //Associative residual chemical potential for CPA
     lnX = X.array().log();
     u_assoc_1 = one_4nc.transpose()*(((lnX.asDiagonal())*one_4c));
     u_assoc_2 = ((dlng_dn)*one_4.transpose()).transpose();
     VectorXd u_assoc_2_vector(Map<VectorXd>(u_assoc_2.data(), u_assoc_2.cols()*u_assoc_2.rows()));
     u_assoc_3 = u_assoc_2_vector.transpose()*(one_4nc-X);
     u_assoc_3b = n_v.array()*u_assoc_3;
     u_assoc_4 = (one_nc.transpose())*u_assoc_3b;
     u_assoc = u_assoc_1.array()-0.5*u_assoc_4;
/*
     cout << "0.5*((one_nc.transpose()*u_assoc_4).array()) = \n" << 0.5*((one_nc.transpose()*u_assoc_4).array()) << endl;
     cout << "u_assoc_3 = \n" << u_assoc_3 << endl;
     cout << "u_assoc_2 = \n" << u_assoc_2 << endl;
     cout << "u_assoc_1.array()-0.5*((one_nc.transpose()*u_assoc_4).array()) = \n" << u_assoc_1.array()-0.5*((one_nc.transpose()*u_assoc_4).array()) << endl;
     cout << "u_assoc_1.array() = \n" << u_assoc_1.array() << endl;
     cout << "u_assoc_2_vector = \n" << u_assoc_2_vector << endl;
     cout << "u_assoc = \n" << u_assoc << endl;;
*/

     //Physical(SRK) residual chemical potential for CPA
     Fn = -log(1-Bcpa/Vt);
     f = (log(1+Bcpa/Vt))/(R*Bcpa);
     f = 1/(R*Bcpa)*(log(1+Bcpa/Vt));
     fV = -1/(R*Vt*(Vt+Bcpa));
     D_T = n_v.transpose()*(aij*n_v);
     FB = n_t/(Vt-Bcpa)+D_T*(f+Vt*fV)/(T*Bcpa);
     FD = -f/T;
     u_SRK = Fn + FB*(Biv.array()) + FD*(Div.array());

     //Cálculo do potencial químico residual da CPA
     u_CPA = u_assoc + u_SRK;
/*
     cout << "V = " << V << endl;

     cout << "Bcpa = " << Bcpa << endl;
     cout << "Biv = " << Biv << endl;
     cout << "Div = " << Div << endl;
     cout << "FB = " << FB << endl;
     cout << "FD = " << FD << endl;
     cout << "D_T = " << D_T << endl;
     cout << "fV = " << fV << endl;
     cout << "f = " << f << endl;
     cout << "Fn = " << Fn << endl;

     cout << "Z_SRK = " << Z_SRK << endl;
     cout << "Z_assoc = " << Z_assoc << endl;
     cout << "Z_CPA = " << Z_CPA << endl;
     cout << "log Z_CPA = " << log(Z_CPA) << endl;
     cout << "u_assoc = " << u_assoc << endl;
     cout << "u_SRK = " << u_SRK << endl;
*/

return u_CPA;
    }

/*
//Função com função objetiva para cálculo de V para CPA
double CPA_volume_obj_function(int nc, double V, double R, double P, double am, double bm, double T, VectorXd x, VectorXd X, double B)
{
    VectorXd one_four(4), one_nc(nc);
    double F;
    int i;

    one_four << 1,
                1,
                1,
                1;

    for(i=0; i<nc; i++)
    {
         one_nc(i) = 1;
    }

    F = R*T/(V-bm) - am/(V*(V+bm)) - 0.5*R*T/V * (1+0.475*B/(V-0.475*B)) * (((x*one_four.transpose()).transpose()).asDiagonal()*one_nc).transpose() * X - P;
    return F;
}

//Função para calcular o volume
double volume_function(int nc, int EdE, int phase, VectorXd x, VectorXd X, VectorXd EdE_parameters, double bm,
                       double am, double R, double T, double P, double tolV, double tolZ, VectorXd b)
{
    int i;
    double B, sigma, epsilon, OMEGA, PSI, Vi, errorV, V, errorZ, condV, Vnew, F_obj, F_obj_plus, F_obj_minus;
    double Zi, Z, qe;

    tolZ = tolV;

B = bm*P/(R*T);
qe = (am/(bm*R*T));



sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    switch(EdE)
    {
    case 1: //SRK
        switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Zi = B; //Chute inicial pra fase líquida
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;


        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Zi = 1; //Chute inicial pra fase vapor
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;
        }

    case 2: //PR
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Zi = B; //Chute inicial pra fase líquida
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;


        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Zi = 1; //Chute inicial pra fase vapor
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;
        }

    case 3: //CPA
    //Chute inicial de V parte de SRK
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Zi = B; //Chute inicial pra fase líquida
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = B + (Zi+epsilon*B) * (Zi+sigma*B) * ((1+B-Zi)/(qe*B));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;

        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Zi = 1; //Chute inicial pra fase vapor
        errorZ = tolZ + 1;
        while(errorZ>tolZ)
        {
            Z = 1+ B - (qe*B)*((Zi-B)/((Zi+epsilon*B)*(Zi+sigma*B)));
            errorZ = fabs(Z-Zi);
            Zi = Z;
        }
        V = R*T*Z/P;
        break;
        }

    //Cálculo usando V de SRK como chute inicial---------
    double F_obj, F_obj_minus, F_obj_plus, F_obj_derivate, Bcpa;
    int i;
    VectorXd one(nc), n_v(nc);
    MatrixXd bij(nc,nc);

    //***********************************************************
    n_v = x; //NESTE CASO!!!!!!!!!!
    //***********************************************************

    for(i=0;i<nc;i++)
    {
        one(i) = 0;
    }

    bij = ((b*one.transpose())+((b*one.transpose()).transpose()))/2;
    Bcpa = n_v.transpose()*(bij.diagonal());

    condV = tolV+1;
    while(condV>tolV)
    {
    F_obj = CPA_volume_obj_function(nc, V, R, P, am, bm, T, x, X, Bcpa);
    cout << "F_obj = " << F_obj << endl;
    double e;
    e = 0.001;
    F_obj_plus = CPA_volume_obj_function(nc, V+e, R, P, am, bm, T, x, X, Bcpa);
    F_obj_minus = CPA_volume_obj_function(nc, V-e, R, P, am, bm, T, x, X, Bcpa);

    F_obj_derivate = (F_obj_plus-F_obj_minus)/(2*e);

    condV = fabs(F_obj/F_obj_derivate);

    Vnew = V - F_obj/F_obj_derivate;
    V = Vnew;
    }
    break;

    }

return V;
}


//**************************************************


/* Função desativada, a função ativa para calcular o volume calcula primeiro o fator de compressibilidade e depois converte para V por V = RTZ/P
//Função para calcular o volume
double volume_function(int nc, int EdE, int phase, VectorXd x, VectorXd X, VectorXd EdE_parameters, double bm,
                       double am, double R, double T, double P, double tolV, double alfa)
{
    int i;
    double B, sigma, epsilon, OMEGA, PSI, Vi, errorV, V, errorZ, condV, Vnew, F_obj, F_obj_plus, F_obj_minus;

B = bm*P/(R*T);

sigma = EdE_parameters[0];
epsilon = EdE_parameters[1];
OMEGA = EdE_parameters[2];
PSI = EdE_parameters[3];

    switch(EdE)
    {
    case 1: //SRK
        switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Vi = bm; //Chute inicial pra fase líquida
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = bm + (Vi+epsilon*bm)*(Vi+sigma*bm)*(R*T+bm*P-Vi*P)/alfa;
            errorV = fabs(V-Vi);
            Vi = V;
            }
        break;


        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Vi = R*T/P; //Chute inicial pra fase vapor
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = R*T/P + bm - (am/P)*(Vi-bm)/((Vi+epsilon*bm)*(Vi+sigma*bm));
            errorZ = fabs(V-Vi);
            Vi = V;
            }
        break;
        }
    break;

    case 2: //PR
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Vi = bm; //Chute inicial pra fase líquida
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = bm + (Vi+epsilon*bm)*(Vi+sigma*bm)*(R*T+bm*P-Vi*P)/alfa;
            errorV = fabs(V-Vi);
            Vi = V;
            }
        break;


        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Vi = R*T/P; //Chute inicial pra fase vapor
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = R*T/P + bm - (am/P)*(Vi-bm)/((Vi+epsilon*bm)*(Vi+sigma*bm));
            errorZ = fabs(V-Vi);
            Vi = V;
            }
        break;
        }
    break;

    case 3: //CPA
    //Chute inicial de V parte de SRK
    switch(phase)
        {
        case 1: //FASE LÍQUIDA --------------------------------------------------------------------------------------
        Vi = bm; //Chute inicial pra fase líquida
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = bm + (Vi)*(Vi+bm)*(R*T+bm*P-Vi*P)/alfa;
            errorV = fabs(V-Vi);
            Vi = V;
            }
        break;

        case 2: //FASE VAPOR ----------------------------------------------------------------------------------------
        Vi = R*T/P; //Chute inicial pra fase vapor
        errorV = tolV + 1;
            while(errorV>tolV)
            {
            V = R*T/P + bm - (am/P)*(Vi-bm)/(Vi*(Vi+bm));
            errorZ = fabs(V-Vi);
            Vi = V;
            }
        break;
        }

    //Cálculo usando V de SRK como chute inicial---------
    double F_obj, F_obj_minus, F_obj_plus, F_obj_derivate;

    condV = tolV+1;
    while(condV>tolV)
    {
    F_obj = CPA_volume_obj_function(nc, V, R, P, am, bm, T, x, X, B);

    double e;
    e = 0.001;
    F_obj_plus = CPA_volume_obj_function(nc, V+e, R, P, am, bm, T, x, X, B);
    F_obj_minus = CPA_volume_obj_function(nc, V-e, R, P, am, bm, T, x, X, B);

    F_obj_derivate = (F_obj_plus-F_obj_minus)/(2*e);

    condV = F_obj/F_obj_derivate;

    Vnew = V - F_obj/F_obj_derivate;
    }
    V = Vnew;
    break;

    }

return V;
}
*/

#endif // EDE_H_INCLUDED
