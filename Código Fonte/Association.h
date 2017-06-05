#ifndef ASSOCIATION_H_INCLUDED
#define ASSOCIATION_H_INCLUDED

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

#include "Gibbs.h"
#include "EdE.h"

using namespace std;
using namespace Eigen;

MatrixXd volume_auto_association(int assoc_scheme, double beta)
{

double beta_AA, beta_AB, beta_BB, beta_BA, beta_AC, beta_BC, beta_CC, beta_CB, beta_CA, beta_AD,
       beta_BD, beta_CD, beta_DD, beta_DC, beta_DB, beta_DA;
MatrixXd beta_auto(4,4);

switch(assoc_scheme)
{
case 1:  //1
//Volume de associação
    beta_AA = beta;
    beta_AB = 0;
    beta_BB = 0;
    beta_BA = 0;
    beta_AC = 0;
    beta_BC = 0;
    beta_CC = 0;
    beta_CB = 0;
    beta_CA = 0;
    beta_AD = 0;
    beta_BD = 0;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = 0;
    beta_DC = 0;
    break;

case 2: //2A
//Volume de associação
    beta_AA = beta;
    beta_AB = beta;
    beta_BB = beta;
    beta_BA = beta;
    beta_AC = 0;
    beta_BC = 0;
    beta_CC = 0;
    beta_CB = 0;
    beta_CA = 0;
    beta_AD = 0;
    beta_BD = 0;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = 0;
    beta_DC = 0;
    break;

case 3: //2B
//Volume de associação
    beta_AA = 0;
    beta_AB = beta;
    beta_BB = 0;
    beta_BA = beta;
    beta_AC = 0;
    beta_BC = 0;
    beta_CC = 0;
    beta_CB = 0;
    beta_CA = 0;
    beta_AD = 0;
    beta_BD = 0;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = 0;
    beta_DC = 0;
    break;

case 4: //3A
//Volume de associação
    beta_AA = beta;
    beta_AB = beta;
    beta_BB = beta;
    beta_BA = beta;
    beta_AC = beta;
    beta_BC = beta;
    beta_CC = beta;
    beta_CB = beta;
    beta_CA = beta;
    beta_AD = 0;
    beta_BD = 0;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = 0;
    beta_DC = 0;
    break;

case 5: //3B
//Volume de associação
    beta_AA = 0;
    beta_AB = 0;
    beta_BB = 0;
    beta_BA = 0;
    beta_AC = beta;
    beta_BC = beta;
    beta_CC = 0;
    beta_CB = beta;
    beta_CA = beta;
    beta_AD = 0;
    beta_BD = 0;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = 0;
    beta_DC = 0;
    break;

case 6: //4A
//Volume de associação
    beta_AA = beta;
    beta_AB = beta;
    beta_BB = beta;
    beta_BA = beta;
    beta_AC = beta;
    beta_BC = beta;
    beta_CC = beta;
    beta_CB = beta;
    beta_CA = beta;
    beta_AD = beta;
    beta_BD = beta;
    beta_CD = beta;
    beta_DD = beta;
    beta_DA = beta;
    beta_DB = beta;
    beta_DC = beta;
    break;

case 7: //4B
//Volume de associação
    beta_AA = 0;
    beta_AB = 0;
    beta_BB = 0;
    beta_BA = 0;
    beta_AC = 0;
    beta_BC = 0;
    beta_CC = 0;
    beta_CB = 0;
    beta_CA = 0;
    beta_AD = beta;
    beta_BD = beta;
    beta_CD = beta;
    beta_DD = 0;
    beta_DA = beta;
    beta_DB = beta;
    beta_DC = beta;
    break;

/*
case 8: //4C
//Volume de associação
    beta_AA = 0;
    beta_AB = 0;
    beta_BB = 0;
    beta_BA = 0;
    beta_AC = beta;
    beta_BC = beta;
    beta_CC = 0;
    beta_CB = beta;
    beta_CA = beta;
    beta_AD = beta;
    beta_BD = beta;
    beta_CD = 0;
    beta_DD = 0;
    beta_DA = beta;
    beta_DB = beta;
    beta_DC = 0;
    break;
*/

case 8: //4C
//Volume de associação
    beta_AA = 0;
    beta_AB = beta;
    beta_BB = 0;
    beta_BA = beta;
    beta_AC = beta;
    beta_BC = 0;
    beta_CC = 0;
    beta_CB = 0;
    beta_CA = beta;
    beta_AD = 0;
    beta_BD = beta;
    beta_CD = beta;
    beta_DD = 0;
    beta_DA = 0;
    beta_DB = beta;
    beta_DC = beta;
    break;
}

//Define a matriz beta para volume de auto-associação
beta_auto << beta_AA, beta_AB, beta_AC, beta_AD,
             beta_BA, beta_BB, beta_BC, beta_BD,
             beta_CA, beta_CB, beta_CC, beta_CD,
             beta_DA, beta_DB, beta_DC, beta_DD;

return beta_auto;
}

//Função para calcular a energia de auto-associação-------------------------------------------------
MatrixXd energy_auto_association(int assoc_scheme, double E, double R, double T)
{
double E_AA, E_AB, E_BB, E_BA, E_AC, E_BC, E_CC, E_CB, E_CA, E_AD, E_BD, E_CD, E_DD, E_DC, E_DB, E_DA;
MatrixXd E_auto(4,4), E_auto_assoc(4,4);

switch(assoc_scheme)
{
case 1:  //1
//Energia de associação
    E_AA = E;
    E_AB = 0;
    E_BB = 0;
    E_BA = 0;
    E_AC = 0;
    E_BC = 0;
    E_CC = 0;
    E_CB = 0;
    E_CA = 0;
    E_AD = 0;
    E_BD = 0;
    E_CD = 0;
    E_DD = 0;
    E_DA = 0;
    E_DB = 0;
    E_DC = 0;
    break;

case 2: //2A
//Energia de associação
    E_AA = E;
    E_AB = E;
    E_BB = E;
    E_BA = E;
    E_AC = 0;
    E_BC = 0;
    E_CC = 0;
    E_CB = 0;
    E_CA = 0;
    E_AD = 0;
    E_BD = 0;
    E_CD = 0;
    E_DD = 0;
    E_DA = 0;
    E_DB = 0;
    E_DC = 0;
    break;

case 3: //2B
//Energia de associação
    E_AA = 0;
    E_AB = E;
    E_BB = 0;
    E_BA = E;
    E_AC = 0;
    E_BC = 0;
    E_CC = 0;
    E_CB = 0;
    E_CA = 0;
    E_AD = 0;
    E_BD = 0;
    E_CD = 0;
    E_DD = 0;
    E_DA = 0;
    E_DB = 0;
    E_DC = 0;
    break;

case 4: //3A
//Energia de associação
    E_AA = E;
    E_AB = E;
    E_BB = E;
    E_BA = E;
    E_AC = E;
    E_BC = E;
    E_CC = E;
    E_CB = E;
    E_CA = E;
    E_AD = 0;
    E_BD = 0;
    E_CD = 0;
    E_DD = 0;
    E_DA = 0;
    E_DB = 0;
    E_DC = 0;
    break;

case 5: //3B
//Energia de associação
    E_AA = 0;
    E_AB = 0;
    E_BB = 0;
    E_BA = 0;
    E_AC = E;
    E_BC = E;
    E_CC = 0;
    E_CB = E;
    E_CA = E;
    E_AD = 0;
    E_BD = 0;
    E_CD = 0;
    E_DD = 0;
    E_DA = 0;
    E_DB = 0;
    E_DC = 0;
    break;

case 6: //4A
//Energia de associação
    E_AA = E;
    E_AB = E;
    E_BB = E;
    E_BA = E;
    E_AC = E;
    E_BC = E;
    E_CC = E;
    E_CB = E;
    E_CA = E;
    E_AD = E;
    E_BD = E;
    E_CD = E;
    E_DD = E;
    E_DA = E;
    E_DB = E;
    E_DC = E;
    break;

case 7: //4B
//Energia de associação
    E_AA = 0;
    E_AB = 0;
    E_BB = 0;
    E_BA = 0;
    E_AC = 0;
    E_BC = 0;
    E_CC = 0;
    E_CB = 0;
    E_CA = 0;
    E_AD = E;
    E_BD = E;
    E_CD = E;
    E_DD = 0;
    E_DA = E;
    E_DB = E;
    E_DC = E;
    break;
/*
case 8: //4C
//Energia de associação
    E_AA = 0;
    E_AB = 0;
    E_BB = 0;
    E_BA = 0;
    E_AC = E;
    E_BC = E;
    E_CC = 0;
    E_CB = E;
    E_CA = E;
    E_AD = E;
    E_BD = E;
    E_CD = 0;
    E_DD = 0;
    E_DA = E;
    E_DB = E;
    E_DC = 0;
    break;
*/
case 8: //4C
//Energia de associação
    E_AA = 0;
    E_AB = E;
    E_BB = 0;
    E_BA = E;
    E_AC = E;
    E_BC = 0;
    E_CC = 0;
    E_CB = 0;
    E_CA = E;
    E_AD = 0;
    E_BD = E;
    E_CD = E;
    E_DD = 0;
    E_DA = 0;
    E_DB = E;
    E_DC = E;
    break;
}

//Define a matriz E para energia de auto-associação
E_auto << E_AA, E_AB, E_AC, E_AD,
          E_BA, E_BB, E_BC, E_BD,
          E_CA, E_CB, E_CC, E_CD,
          E_DA, E_DB, E_DC, E_DD;

E_auto_assoc = (E_auto.array()/(R*T)).exp();

return E_auto_assoc;
}

//Função para calcular o volume de auto-associação--------------------------------------------------
MatrixXd volume_cross_association(int assoc_scheme, double beta)
{

double beta_AA, beta_AB, beta_BB, beta_BA, beta_AC, beta_BC, beta_CC, beta_CB, beta_CA, beta_AD,
       beta_BD, beta_CD, beta_DD, beta_DC, beta_DB, beta_DA;
MatrixXd beta_auto(4,2);

switch(assoc_scheme)
{
case 1:  //1
//Volume de associação
    beta_auto << beta, beta,
                 0,    0,
                 0,    0,
                 0,    0;
    break;

case 2: //2A
//Volume de associação
    beta_auto << beta,    beta,
                 beta,    beta,
                 0,    0,
                 0,    0;
    break;

case 3: //2B
//Volume de associação
    beta_auto << 0,    beta,
                 beta, 0,
                 0,    0,
                 0,    0;
    break;

case 4: //3A
//Volume de associação
    beta_auto << beta,    beta,
                 beta,    beta,
                 beta,    beta,
                 0,    0;
    break;

case 5: //3B
//Volume de associação
    beta_auto << 0,    beta,
                 0,    beta,
                 beta,    0,
                 0,       0;
    break;

case 6: //4A
//Volume de associação
    beta_auto << beta,    beta,
                 beta,    beta,
                 beta,    beta,
                 beta,    beta;
    break;

case 7: //4B
//Volume de associação
    beta_auto << 0,    0,
                 0,    0,
                 0,    0,
                 0,    0;
    break;


case 8: //4C
//Volume de associação
    beta_auto << 0,    beta,
                 0,    beta,
                 beta,    0,
                 beta,    0;
    break;
}

return beta_auto;
}

//Função para calcular a energia de auto-associação-------------------------------------------------
MatrixXd energy_cross_association(int assoc_scheme, double E, double R, double T)
{
double E_AA, E_AB, E_BB, E_BA, E_AC, E_BC, E_CC, E_CB, E_CA, E_AD, E_BD, E_CD, E_DD, E_DC, E_DB, E_DA;
MatrixXd E_auto(4,2), E_auto_assoc(4,2);

switch(assoc_scheme)
{
case 1:  //1
//Energia de associação
    E_auto << E, E,
              0, 0,
              0, 0,
              0, 0;
    break;

case 2: //2A
//Energia de associação
    E_auto << 0, 0,
              0, 0,
              0, 0,
              0, 0;
    break;

case 3: //2B
//Energia de associação
    E_auto << 0, E,
              E, 0,
              0, 0,
              0, 0;
    break;

case 4: //3A
//Energia de associação
    E_auto << 0, 0,
              0, 0,
              0, 0,
              0, 0;
    break;

case 5: //3B
//Energia de associação
    E_auto << 0, E,
              0, E,
              E, 0,
              0, 0;
    break;

case 6: //4A
//Energia de associação
    E_auto << E, E,
              E, E,
              E, E,
              E, E;
    break;

case 7: //4B
//Energia de associação
    E_auto << 0, 0,
              0, 0,
              0, 0,
              0, 0;
    break;

case 8: //4C
//Energia de associação
    E_auto << 0, E,
              0, E,
              E, 0,
              E, 0;
    break;
}

E_auto_assoc = (E_auto.array()/(2*R*T)).exp();

return E_auto_assoc;
}

//Função para calcular a matriz DELTA---------------------------------------------------------------
MatrixXd DELTA_function(int combining_rule, int nc, int phase, double R, double T, double P, double tolV,
               VectorXd alfa, double am, double bm, MatrixXd beta_col, MatrixXd beta_row, MatrixXd E_col,
               MatrixXd E_row, int EdE, VectorXd x, VectorXd X, VectorXd EdE_parameters, VectorXd bi, double tolZ,
               double V, double BETCR, MatrixXd E_auto, MatrixXd beta_auto, double Dij)
{
    int nc4;
    nc4 = 4*nc;

    MatrixXd bij(nc,nc), Bij(nc4,nc4);
    MatrixXd bibj(nc,nc), BiBj2(nc,nc), BiBj(nc4,nc4);

//Calculando a densidade reduzida
double eta;
eta = bm/(4*V);

//Calculando a função de distribuição radial (para sCPA)
double g_Vm;
g_Vm = 1/(1-1.9*eta);

int i;
VectorXd one(nc);
MatrixXd one_4(4,4);

for(i=0;i<nc;i++)
{
    one(i) = 1;
}

one_4 << 1,1,1,1,
         1,1,1,1,
         1,1,1,1,
         1,1,1,1;

bij = ((bi*(one.transpose()))+((bi*(one.transpose())).transpose()))/2;
Bij = kroneckerProduct(bij,one_4);

bibj = ((bi*(one.transpose())).cwiseProduct((bi*(one.transpose())).transpose()));
BiBj = kroneckerProduct(bibj,one_4);
BiBj2 = BiBj.array().pow(0.5);

/*
cout << "b = \n" << bi << endl;
cout << "bij = \n" << bij << endl;
cout << "Bij = \n" << Bij << endl;
cout << "bibj = \n" << bibj << endl;
cout << "BiBj = \n" << BiBj << endl;
cout << "BiBj2 = \n" << BiBj2 << endl;
*/

//Definindo a matriz para energia e volume de associação-cruzada
MatrixXd E_cross(nc4,nc4), E_cross_1(nc4,nc4), beta_cross(nc4,nc4), DELTA(nc4,nc4), DELTA_row(nc4,nc4),
         E_row_1(nc4,nc4), E_col_1(nc4,nc4), DELTA_col(nc4,nc4), beta_cross_BETCR(nc4,nc4), zero4(nc4,nc4),
         E_auto_1(nc4,nc4), DELTA_auto(nc4,nc4), DELTA_cross(nc4,nc4), D12(nc4,nc4);
         int j;

  switch(combining_rule)
  {
  case 1: //CR-1
    zero4 = MatrixXd::Zero(4,4);

    for(i=0;i<4*nc;i++)
    {
        for(j=0;j<2;j++)
        {
            if(E_row(i,j)==1)
            {
                E_row(i,j) = 0;
            }
        }
    }

    E_cross = E_row*E_row.transpose();
    beta_cross = (beta_row*beta_row.transpose()).array().pow(0.5);
/*
    //beta_cross(0,5) = pow((beta_cross(0,1)*beta_cross(4,4)),(0.5));
    beta_cross(1,4) = pow((beta_cross(0,1)*beta_cross(4,4)),(0.5));
    //beta_cross(5,0) = pow((beta_cross(0,1)*beta_cross(4,4)),(0.5));
    beta_cross(4,1) = pow((beta_cross(0,1)*beta_cross(4,4)),(0.5));

    E_cross(1,4) = pow(E_cross(0,1)*E_cross(4,4),0.5);
    E_cross(4,1) = pow(E_cross(0,1)*E_cross(4,4),0.5);
*/
    for(j=0;j<nc;j++)
    {
        E_cross.block(4*j, 4*j, 4, 4) << zero4;
        beta_cross.block(4*j, 4*j, 4, 4) = zero4;
    }

    E_cross = E_cross.array() + E_auto.array();
    beta_cross = beta_cross.array() + beta_auto.array();

    for(i=0;i<4*nc;i++)
    {
        for(j=0;j<4*nc;j++)
        {
            if(E_cross(i,j)<1e-100)
            {
                E_cross(i,j) = 1;
            }
        }
    }

    E_cross_1 = E_cross.array()-1;
    DELTA = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
/*
    cout << "\nE_row = \n" << E_row << endl;
    cout << "\nE_auto = \n" << E_auto << endl;
    cout << "\nbeta_row = \n" << beta_row << endl;
    cout << "\nbeta_auto = \n" << beta_auto << endl;
    cout << "\nE_cross = \n" << E_cross << endl;
    cout << "\nbeta_cross = \n" << beta_cross << endl;
    cout << "\nDELTA = \n" << DELTA << endl;
    cin.get();
*/
    break;

  case 2: //CR-2
    E_cross = (E_row + E_col)/2;
    beta_cross = (beta_row + beta_col)/2;
    E_cross_1 = E_cross.array()-1;
    DELTA = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
    break;

  case 3: //CR-3
    E_cross = (E_row.cwiseProduct(E_col)).array().pow(0.5);
    beta_cross = (beta_row.cwiseProduct(beta_col)).array().pow(0.5);
    E_cross_1 = E_cross.array()-1;
    DELTA = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
    break;

  case 4: //CR-4
    E_cross = (E_row.cwiseProduct(E_col)).array().pow(0.5);
    beta_cross = (beta_row + beta_col)/2;;
    E_cross_1 = E_cross.array()-1;
    DELTA = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
    break;


  case 5: //ECR
    zero4 = MatrixXd::Zero(4,4);

    for(i=0;i<4*nc;i++)
    {
        for(j=0;j<2;j++)
        {
            if(E_row(i,j)==1)
            {
                E_row(i,j) = 0;
            }
        }
    }

    E_cross = E_row*E_row.transpose();
    beta_cross = (beta_row*beta_row.transpose()).array().pow(0.5);

    for(j=0;j<nc;j++)
    {
        E_cross.block(4*j, 4*j, 4, 4) << zero4;
        beta_cross.block(4*j, 4*j, 4, 4) = zero4;
    }

    for(i=0;i<4*nc;i++)
    {
        for(j=0;j<4*nc;j++)
        {
            if(E_cross(i,j)<1e-100)
            {
                E_cross(i,j) = 1;
            }
        }
    }

    E_cross_1 = E_cross.array()-1;
    E_auto_1 = E_auto.array()-1;

    DELTA_auto = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
    DELTA_cross = (g_Vm * E_auto_1).cwiseProduct(Bij.cwiseProduct(beta_auto));
    DELTA = DELTA_auto.array() + DELTA_cross.array();


    break;

    case 6: //mCR-1
    E_cross = E_row.cwiseProduct(E_col);
    beta_cross = (beta_row.cwiseProduct(beta_col)).array().pow(0.5);

    beta_cross_BETCR << 0.00, 0.00, 0.00, 0.00, 0.00, BETCR, 0.00, 0.00,
                        0.00, 0.00, 0.00, 0.00, BETCR, 0.00, 0.00, 0.00,
                        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                        0.00, BETCR, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                        BETCR, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                        0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00;

    //cout << "beta_cross = \n" << beta_cross << endl;
    //cout << "beta_cross_BETCR = \n" << beta_cross_BETCR << endl;

    beta_cross = beta_cross.array()+beta_cross_BETCR.array();

    //cout << "beta_cross final = \n" << beta_cross << endl;

    E_cross_1 = E_cross.array()-1;
    DELTA = (g_Vm * E_cross_1).cwiseProduct(Bij.cwiseProduct(beta_cross));
    break;

  }

  D12 << 1.00, 1.00, 1.00, 1.00, Dij, Dij, Dij, Dij,
         1.00, 1.00, 1.00, 1.00, Dij, Dij, Dij, Dij,
         1.00, 1.00, 1.00, 1.00, Dij, Dij, Dij, Dij,
         1.00, 1.00, 1.00, 1.00, Dij, Dij, Dij, Dij,
         Dij, Dij, Dij, Dij, 1.00, 1.00, 1.00, 1.00,
         Dij, Dij, Dij, Dij, 1.00, 1.00, 1.00, 1.00,
         Dij, Dij, Dij, Dij, 1.00, 1.00, 1.00, 1.00,
         Dij, Dij, Dij, Dij, 1.00, 1.00, 1.00, 1.00;

DELTA = DELTA.cwiseProduct(D12);
return DELTA;
}

//Função INICIAL para calcular a fração de moléculas 'i' com sítios A,B,C e D não-ligados
VectorXd fraction_nbs_initial(int nc, int combining_rule, int phase, double R, double T, double P, double tolV,
                      VectorXd alfa, double am, double bm, MatrixXd beta_col, MatrixXd beta_row, MatrixXd E_col,
                      MatrixXd E_row, double tolX, VectorXd x, int EdE, VectorXd EdE_parameters, VectorXd bi,
                      double tolZ, double V, double BETCR, MatrixXd E_auto, MatrixXd beta_auto, double Dij)
{
    VectorXd X(4*nc), Xnew(4*nc), Xcond(4*nc), pre_Xnew(4*nc), g(4*nc), deltaX(4*nc), lnXX(4*nc), X1(4*nc), lnXX2(4*nc);
    MatrixXd DELTA(4*nc,4*nc), xXD (4*nc,4*nc), one_4_x(4,nc), H(4*nc,4*nc), K(4*nc,4*nc), I(4*nc,4*nc), H_1(4*nc,4*nc);
    MatrixXd X_invmat(4*nc,4*nc), zero_one(4*nc,4*nc), DELTA1(4*nc,4*nc);
    double X_max, Q, Qold;
    int i;

    //Chute inicial de X, todos os sítios com a mesma probabilidade de estarem ocupados
    for(i=0;i<(4*nc);i++)
    {
        X(i) = 0.2;
    }

    VectorXd one(nc);
    for(i=0;i<nc;i++)
    {
        one(i) = 1;
    }

    VectorXd one_4nc(4*nc);
    for(i=0;i<4*nc;i++)
    {
        one_4nc(i) = 1;
    }

    VectorXd one_4(4);
    one_4 << 1, 1, 1, 1;

/*
X_max = tolX + 1;
tolX = 0.001;

cout << "before X loop" << endl;
while(X_max>tolX)
{
    cout << "inside X loop" << endl;
    V = volume_function(nc, EdE, phase, x, X, EdE_parameters, bm, am, R, T, P, tolV, tolZ);
    cout << "V = " << V << endl;
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ);
    cout << "DELTA = \n" << DELTA << endl;
    //xXD = (((x.asDiagonal()*X)*one_4nc).transpose()).cwiseProduct(DELTA);
    //Xnew = one_4nc + 1/V*(one_4nc.transpose()*xXD);
    cout << "x.asDiagonal()*X = \n" << x.asDiagonal()*X << endl;
    cout << "xXD = \n" << xXD << endl;
    cout << "Xnew = \n" << Xnew << endl;
    cout << "X = \n" << X << endl;

    Xcond = Xnew - X;
    cout << "Xcond = " << Xcond << endl;
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    X = Xnew;
}
*/

one_4_x = one_4*x.transpose();
    VectorXd one_4_x_vector(Map<VectorXd>(one_4_x.data(), one_4_x.cols()*one_4_x.rows()));

I.setIdentity();
i = 1;
tolX = 0.000001;
X_max = tolX + 1;

/*
while(X_max>tolX)
{
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ, V);
    one_4_x = one_4*x.transpose();
    VectorXd one_4_x_vector(Map<VectorXd>(one_4_x.data(), one_4_x.cols()*one_4_x.rows()));
    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/V;
    Xnew = one_4nc + pre_Xnew;
    Xnew = Xnew.asDiagonal().inverse().diagonal();
    //Xnew = 0.8*Xnew+0.2*X;

    Xcond = Xnew - X;
    Xcond = Xcond.array().abs();
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    X = Xnew;
}
*/
while(X_max>tolX)
{
    if(i<=5) //Os primeiros cinco passos são dados com o método de substituição sucessiva
    {
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ, V, BETCR, E_auto, beta_auto, Dij);

    zero_one << 0, 0, 0, 0, 1, 1, 1, 1,
               0, 0, 0, 0, 1, 1, 1, 1,
               0, 0, 0, 0, 1, 1, 1, 1,
               0, 0, 0, 0, 1, 1, 1, 1,
               1, 1, 1, 1, 0, 0, 0, 0,
               1, 1, 1, 1, 0, 0, 0, 0,
               1, 1, 1, 1, 0, 0, 0, 0,
               1, 1, 1, 1, 0, 0, 0, 0;
    DELTA1 = DELTA.cwiseProduct(zero_one);


    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/V;
    Xnew = one_4nc + pre_Xnew;
    Xnew = Xnew.asDiagonal().inverse().diagonal();
    Xnew = 0.8*Xnew+0.2*X;

    Xcond = Xnew - X;
    Xcond = Xcond.array().abs();
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    X = Xnew;

    }

    else //Os passos restantes são dados com o método de segunda ordem
    {
    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/V;
    H_1 = ((X.asDiagonal().inverse())*(one_4_x_vector+(X.transpose()*K).transpose())).asDiagonal();
    //*************
    //H_1 = ((X.asDiagonal().inverse())*(one_4_x_vector+(K*X))).asDiagonal();
    H = -H_1-K;
    X_invmat = X.asDiagonal().inverse();
    g = ((X_invmat-I)*one_4_x_vector).transpose()-X.transpose()*K;
    //*************
    //g = ((X_invmat-I)*one_4_x_vector).transpose()-K*X;

    deltaX = H.inverse()*(-g);

    lnXX = X.array().log()-X.array();
    X1 = lnXX.array()+1;

    Qold =  X.transpose()*((X.transpose()*K.transpose()).transpose());
    Qold = one_4_x_vector.transpose()*(X1)-0.5*Qold;
    //************
    //Qold =  X.transpose()*((X.transpose()*K.transpose()).transpose());
    //Qold = one_4_x_vector.transpose()*(X1)-0.5*Qold;
    Q = Qold-1; //Valor arbitrário para forçar entrada no loop da função objetivo

    while(Q<Qold)
    {
    Xnew = X + deltaX;

    int l;
    for(l=0;l<4*nc;l++)
    {
    Xnew(l) = max(Xnew(l),0.2*(X(l)));
    }
    lnXX = Xnew.array().log()-Xnew.array();
    X1 = lnXX.array()+1;
    Q =  Xnew.transpose()*((Xnew.transpose()*K.transpose()).transpose());
    Q = one_4_x_vector.transpose()*(X1)-0.5*Q;
    //****************
    //Q =  Xnew.transpose()*((Xnew.transpose()*K.transpose()).transpose());
    //Q = one_4_x_vector.transpose()*(X1)-0.5*Q;

    if(Q<Qold)
    {
        deltaX = deltaX/2;
    }
    }

    Xcond = Xnew - X;
    Xcond = Xcond.array().abs();
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    X = Xnew;
    }
    i++;
}

    return X;
}

//Função para calcular a fração de moléculas 'i' com sítios A,B,C e D não-ligados
VectorXd fraction_nbs(int nc, int combining_rule, int phase, double R, double T, double P, double tolV,
                      VectorXd alfa, double am, double bm, MatrixXd beta_col, MatrixXd beta_row, MatrixXd E_col,
                      MatrixXd E_row, double tolX, VectorXd x, int EdE, VectorXd EdE_parameters,
                      VectorXd bi, double tolZ, double V, double deltaV, VectorXd X, int iter, VectorXd a,
                      double *Q_func, double BETCR, MatrixXd E_auto, MatrixXd beta_auto, double Dij)
{
    VectorXd Xnew(4*nc), Xnew2(4*nc), Xcond(4*nc), pre_Xnew(4*nc), pre_XnewV2(4*nc), g(4*nc), deltaX(4*nc), lnXX(4*nc),
             X1(4*nc), lnXX2(4*nc), dX_dV(4*nc), QXV(4*nc), X2(4*nc);
    MatrixXd DELTA(4*nc,4*nc), xXD (4*nc,4*nc), one_4_x(4,nc), H(4*nc,4*nc), K(4*nc,4*nc), I(4*nc,4*nc), H_1(4*nc,4*nc), QXX(4*nc,4*nc),
             QXX1(4*nc,4*nc);
    MatrixXd X_invmat(4*nc,4*nc), DELTA1(4*nc,4*nc), zero_one(4*nc,4*nc);
    double X_max, Q, Qold;
    int i;

    double DELTA_pure;

    VectorXd one(nc);
    for(i=0;i<nc;i++)
    {
        one(i) = 1;
    }

    VectorXd one_4nc(4*nc);
    for(i=0;i<4*nc;i++)
    {
        one_4nc(i) = 1;
    }

    VectorXd one_4(4);
    one_4 << 1, 1, 1, 1;

    one_4_x = one_4*x.transpose();
    VectorXd one_4_x_vector(Map<VectorXd>(one_4_x.data(), one_4_x.cols()*one_4_x.rows()));

    //Chute inicial de X a partir de deltaV
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ, V, BETCR, E_auto, beta_auto, Dij);

    DELTA_pure = (1/(1-1.9*(bm/4/V)))*(exp(E_col(1,0)/R/T)-1)*((bi(0)+bi(1))/2)*beta_col(0,1);

    zero_one << 0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0;

    zero_one << 1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1;

    DELTA1 = DELTA.cwiseProduct(zero_one);


    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    pre_Xnew = ((one_4nc.transpose()*xXD).array())/V;

    //pre_Xnew = ((one_4nc.transpose()*xXD).array())/(2*V);

    QXV = (one_4_x_vector.asDiagonal()*pre_Xnew)*(1/V-1);
    X2 = X.array().pow(2);
    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/V;
    QXX1 = (((X2.asDiagonal().inverse())*one_4_x_vector).asDiagonal());
    QXX = -QXX1-K;
    dX_dV = QXX.inverse()*(-QXV);
    X = X + dX_dV*deltaV;

    if(iter==0)
    {
    for(i=0;i<(4*nc);i++)
    {
        X(i) = 0.2;
    }
    }

I.setIdentity();
i = 1;
tolX = 0.000001;
X_max = tolX + 1;

while(X_max>tolX || g.maxCoeff()>tolX)
{
    if(i<=4) //Os primeiros cinco passos são dados com o método de substituição sucessiva
    {
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ, V, BETCR, E_auto, beta_auto, Dij);

    zero_one << 0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0;

    zero_one << 1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1;

    DELTA1 = DELTA.cwiseProduct(zero_one);


    xXD = (((one_4_x_vector.asDiagonal()*X)*one_4nc.transpose()).transpose()).cwiseProduct(DELTA);
    //********O TERMO DE pre_Xnew ESTÁ MODIFICADO, ORIGINALMENTE NÃO SE TEM O 2 MULTIPLICANDO V!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    pre_Xnew = ((one_4nc.transpose()*xXD).array())/V;

    //pre_Xnew = ((one_4nc.transpose()*xXD).array())/(2*V);

    Xnew = one_4nc + pre_Xnew;
    Xnew = Xnew.asDiagonal().inverse().diagonal();
    Xnew = 0.8*Xnew+0.2*X;

    Xcond = Xnew - X;
    Xcond = Xcond.array().abs();
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    X_max = X_max/(X.maxCoeff());

    X = Xnew;
    }

    else //Os passos restantes são dados com o método de segunda ordem
    {
    DELTA = DELTA_function(combining_rule, nc, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                           EdE, x, X, EdE_parameters, bi, tolZ, V, BETCR, E_auto, beta_auto, Dij);

    zero_one << 1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 0, 1, 1, 1, 1;

    DELTA1 = DELTA.cwiseProduct(zero_one);

    K = ((one_4_x_vector*one_4_x_vector.transpose()).cwiseProduct(DELTA))/V;
    H_1 = ((X.asDiagonal().inverse())*(one_4_x_vector+(X.transpose()*K).transpose())).asDiagonal();
    //*************
    //H_1 = ((X.asDiagonal().inverse())*(one_4_x_vector+(K*X))).asDiagonal();
    H = -H_1-K;
    X_invmat = X.asDiagonal().inverse();
    g = ((X_invmat-I)*one_4_x_vector).transpose()-X.transpose()*K;
    //*************
    //g = ((X_invmat-I)*one_4_x_vector).transpose()-K*X;

    deltaX = H.inverse()*(-g);

    lnXX = X.array().log()-X.array();
    X1 = lnXX.array()+1;

    Qold =  X.transpose()*((X.transpose()*K.transpose()).transpose());
    Qold = one_4_x_vector.transpose()*(X1)-0.5*Qold;
    //************
    //Qold =  X.transpose()*((X.transpose()*K.transpose()).transpose());
    //Qold = one_4_x_vector.transpose()*(X1)-0.5*Qold;
    Q = Qold-1; //Valor arbitrário para forçar entrada no loop da função objetivo

    while(Q<Qold)
    {
    Xnew = X + deltaX;

    int l;
    for(l=0;l<4*nc;l++)
    {
    Xnew(l) = max(Xnew(l),0.2*(X(l)));
    }
    lnXX = Xnew.array().log()-Xnew.array();
    X1 = lnXX.array()+1;
    Q =  Xnew.transpose()*((Xnew.transpose()*K.transpose()).transpose());
    Q = one_4_x_vector.transpose()*(X1)-0.5*Q;
    //****************
    //Q =  Xnew.transpose()*((Xnew.transpose()*K.transpose()).transpose());
    //Q = one_4_x_vector.transpose()*(X1)-0.5*Q;

    if(Q<Qold)
    {
        deltaX = deltaX/2;
    }
    }

    (*Q_func) = Q;

    Xcond = Xnew - X;
    Xcond = Xcond.array().abs();
    X_max = Xcond.maxCoeff(); //Valor necessário para convergência

    //X_max = X_max/(X.maxCoeff());

    X = Xnew;

    //ESSE IF ESTAVA ANTES DO Q_FUNC
    }

    i++;
}
/*
    X(0) = ((-1+pow((1+4/V*DELTA(0,1)),0.5))/(2/V*DELTA(0,1)));
    X(1) = ((-1+pow((1+4/V*DELTA(0,1)),0.5))/(2/V*DELTA(0,1)));
    X(2) = 1;
    X(3) = 1;
    X(4) = 1;
    X(5) = 1;
    X(6) = 1;
    X(7) = 1;
*/
    //cout << "X_max = " << X_max << endl;

    return X;
}



#endif // ASSOCIATION_H_INCLUDED
