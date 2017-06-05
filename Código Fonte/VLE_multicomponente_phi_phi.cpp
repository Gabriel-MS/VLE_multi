/*
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                ||PROGRAMA PARA CÁLCULO DE EQUILÍBRIO VLE EM MISTURAS MULTICOMPONENTES||
                ||AUTOR: GABRIEL MORAES SILVA                                         ||
                ||LINGUAGEM: C++                                                      ||
                ||BIBLIOTECAS: EIGEN                                                  ||
                ||ANO: 2016                                                           ||
                ||VERSÃO 1.2                                                          ||
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/


#include <cmath>
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>

typedef std::numeric_limits< double > dbl;
typedef std::numeric_limits< long double > ldbl;

#include "mixingrules.h"
#include "EdE.h"
#include "MSA.h"
#include "Gibbs.h"
#include "Association.h"
#include "Renormalization.h"
#include "interpolation_util.h"
#include "numerical.h"
#include "envelope.h"

template<class Function>
double deriv1(const Function& f, double x)
{
    double dx=1e-8*(1.0+fabs(x));
    double x0=x-dx;
    double x2=x+dx;
    return (f(x2)-f(x0)) / (2.0*dx);
}

//define e include necessários para trabalhar com Eigen
#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen; //Biblioteca para trabalhar com matrizes e vetores


int main()
{

//double aaa;
//aaa = fugacity_renormalized(1,0.5,1,0.0000651244168963985,8.314e-6,273.15);



//Arquivo de saída dos dados--------------------------------------------------------------------------
ofstream output("../Planilhas de análise/Output.csv");
output << "Dados da simulação \n -------------------------------------------------------------------------------------" << endl;


//Apresentação e versão do programa
cout << "    |===========================||" << endl;
cout << "    |Autor: Gabriel Moraes Silva||" << endl;
cout << "    |Ano: 2017                  ||" << endl;
cout << "    |V. 3.0                     ||" << endl;
cout << "    |===========================||" << endl;

int nc, i, n, row, col, phase;
int EdE, MR, process, binary_interaction, G_ex_model;
int mixture;

cout.precision(10);

//Usuário escolhe o número de componentes
cout << "\nNumber of components in mixture: ";
cin >> mixture;

switch(mixture)
{
    case 1: nc = 2; break;
    case 2: nc = 2; break;
    case 3:
        cout << "Define number of components: ";
        cin >> nc;
        break;
}
output << "Number of Components = " << nc << endl;

//Usuário escolhe os componentes
int cp[nc];
for(i=0; i<nc; i++)
{
cout << "Choose component " << i+1 << ": ";
cin >> cp[i];
output << "Component " << i+1 << " = " << cp[i] << endl;
}

//Variáveis------------------------------------------------------------------------
VectorXd Tc(nc), Pc(nc), omega(nc), r(nc), q(nc), q_prime(nc), q_(nc), A(nc), B(nc), C(nc), CT(nc), Psat(nc);
VectorXd L_rg(nc), phi_rg(nc); //Renormalization
VectorXd Tr(nc), alfa(nc), a(nc), b(nc), EdE_parameters(nc), K(nc), Kx(nc), one(nc);
VectorXd omega_virtual(nc), Tc_virtual(nc), Pc_virtual(nc);
VectorXd ln10(nc), logPsat(nc), lnPsat(nc), x(nc), y(nc), phi_liquid_phase(nc), phi_vapor_phase(nc);
VectorXd a0(nc), bCPA(nc), E(nc), beta(nc), c1(nc), X(4*nc), Xl(4*nc), Xv(4*nc);
VectorXd n_v(nc), n_vx(nc), n_vy(nc), yinit(nc);
MatrixXd DELTA (4*nc,4*nc), alfa_NRTL(nc,nc), Aij(nc,nc);
MatrixXd E_row(4*nc,2), E_col(4*nc,2), beta_row(4*nc,2), beta_col(4*nc,2), E_i(4,2), beta_i(4,2), E_iRT(4,2);
MatrixXd E_n(4*nc,4*nc), beta_n(4*nc,4*nc), beta_auto(4*nc,4*nc), E_auto(4*nc,4*nc);
VectorXd u_liquid(nc), u_vapor(nc), cond_u(nc);
double tol_u, max_cond_u, u_liquid1, u_vapor1;
int combining_rule, assoc_scheme;
double Tc_a[nc], Pc_a[nc], omega_a[nc], r_a[nc], q_a[nc], q__a[nc], A_a[nc], B_a[nc], C_a[nc]; //Vetores em C++ apenas para gravar valores do arquivo
double L_a[nc], phi_a[nc]; //Renormalization
double E_a[nc], beta_a[nc], a0_a[nc], bCPA_a[nc], c1_a[nc];
double Tc_v[nc], Pc_v [nc], omega_v[nc];
double P, T, R, Pnew, sumKx, al, bl, av, bv, Ey;
double tolZv, tolZl, tolSUMKx, tolKx, initialSUMKx, sumKxnew, finalSUMKx, tolX, tolV;
double errorZv, errorZl, errorSUMKx, errorKx, k12, V, dP_dV, rho_l, X1, Vl, Vv, Vt, deltaV;
double Pinit, Vlinit, Vvinit, Zl, Zv, X1l, X1v;
VectorXd ln_Ki(nc), pre_P1(nc), Tr_1(nc), pre_P1_exp(nc), pre_ln_Ki;
double Vl_obj, Vv_obj, Ql, Qv, dP_dVl, dP_dVv;
double log10P, Tb, Tinit, Told, Dij;
double G_ex, Vl1, Vl2, Vv1, Vv2;
VectorXd Tsat(nc), Alog10P(nc), gama(nc), ln_gama(nc);
double init_T, final_T, step, BETCR, init_P, final_P, Pold;
int max_num_iter, counter, stop, Renormalization, sr_type, r_type, Iteration;

std::vector<double> rho_rv(1000), x_rv(200);
double** d2P;
double** d2P1;
double** d2P2;
double** d2u;
double **Pmat;
double **P1mat;
double **P2mat;
double **umat;
double **u1mat;
double **u2mat;
double **p2;
double **d2u1;
double **d2u2;
double **Wmat;
double **Dmat;

//--------------------------------------------------------------------------------
max_num_iter = 500;

for(i=0; i<nc; i++)
{
//cout << "Moles of component " << i+1 << ": ";
//cin >> n_v[i];
n_v[i] = 1;
}

//DATA INPUT FROM FILES
//Reading Data Bank------------------------
double prop[150][30]; //Matrix to hold all data from properties.csv organized

    ifstream file("../Planilhas de análise/properties.csv");

    for(int row = 0; row < 150; ++row)
    {
        string line;
        getline(file, line);
        if ( !file.good() )
            break;
        stringstream iss(line);

        for (int col = 0; col < 30; ++col)
        {
            string val;
            getline(iss, val, ';');
            if ( !iss )
                break;
            stringstream convertor(val);
            convertor >> prop[row][col];
        }
    }

//Choosing EoS
cout << "\nChoose the Model: \n 1.Soave-Redlich-Kwong \n 2.Peng-Robinson \n 3.CPA-SRK \n 4.MSA" << endl;
cin >> EdE;
output << "Equation of state = " << EdE << endl;

cout << "\nConsider Renormalization? \n 1.Yes \n 2.No \n 3.Use renormalized data from previous simulation" << endl;
cin >> Renormalization;

if(Renormalization==3)
{
    cout << "Renormalization steps taken: " << endl;
    int rstep;
    cin >> rstep;

    rho_rv.resize(rstep);

        d2P = new double *[200];
        for(int k = 0; k <200; k++)
            d2P[k] = new double[rstep];

        d2P1 = new double *[200];
        for(int k = 0; k <200; k++)
            d2P1[k] = new double[rstep];

        d2P2 = new double *[200];
        for(int k = 0; k <200; k++)
            d2P2[k] = new double[rstep];

        d2u1 = new double *[200];
        for(int k = 0; k <200; k++)
            d2u1[k] = new double[rstep];

        d2u2 = new double *[200];
        for(int k = 0; k <200; k++)
            d2u2[k] = new double[rstep];

        Pmat = new double *[200];
        for(int k = 0; k <200; k++)
            Pmat[k] = new double[rstep];

        P1mat = new double *[200];
        for(int k = 0; k <200; k++)
            P1mat[k] = new double[rstep];

        P2mat = new double *[200];
        for(int k = 0; k <200; k++)
            P2mat[k] = new double[rstep];

        umat = new double *[200];
        for(int k = 0; k <200; k++)
            umat[k] = new double[rstep];

        u1mat = new double *[200];
        for(int k = 0; k <200; k++)
            u1mat[k] = new double[rstep];

        u2mat = new double *[200];
        for(int k = 0; k <200; k++)
            u2mat[k] = new double[rstep];


    d2Pgen(d2P,rstep);
    d2P1gen(d2P1,rstep);
    d2P2gen(d2P2,rstep);
    d2u1gen(d2u1,rstep);
    d2u2gen(d2u2,rstep);
    //new_d2u1gen(d2u1,rstep);
    //new_d2u2gen(d2u2,rstep);
    renorm_mat_reader(Pmat,umat,rstep);
    renorm_pp_reader(P1mat,P2mat,rstep);
    renorm_uu_reader(u1mat,u2mat,rstep);
    rho_rv = renorm_rhovec(rstep);
    x_rv = renorm_xvec(rstep);
    d2ugen(d2u1,d2u2,u1mat,u2mat,rstep);
    if(EdE==1) EdE = 5;
    if(EdE==3) EdE = 6;
    Renormalization = 2;
}



if(Renormalization==1)
{
    cout << "short-range type: \n 1. n \n 2. 2n+1 \n 3. 2n-1 \n 4. phi/L^2 \n 5. phi/L^2 (m)" << endl;
    cin >> sr_type;
    cout << "Calculation type: \n 1. Vector \n 2. Point-wise+Vec \n 3. Point-wise" << endl;
    cin >> r_type;
    cout << "Iteration method: \n 1. Trapezoidal (no-function) \n 2. Trapezoidal \n 3. Simpson" << endl;
    cin >> Iteration;
}

output << "Renormalization = " << Renormalization << endl;

//Transferring values from 'prop' matrix to vectors
for(n=0; n<nc; n++)
{
row = cp[n] - 1;
Tc_a[n] = prop[row][2];
Pc_a[n] = prop[row][3];
omega_a[n] = prop[row][4];
r_a[n] = prop[row][5];
q_a[n] = prop[row][6];
q__a[n] = prop[row][7];
A_a[n] = prop[row][8];
B_a[n] = prop[row][9];
C_a[n] = prop[row][10];
a0_a[n] = prop[row][13];
bCPA_a[n] = prop[row][14];
c1_a[n] = prop[row][15];
E_a[n] = prop[row][16];
beta_a[n] = prop[row][17];
Tc_v[n] = prop[row][18];
Pc_v[n] = prop[row][19];
omega_v[n] = prop[row][20];
L_a[n] = prop[row][21];
phi_a[n] = prop[row][22];
}

//Reading C++ vectors into Eigen type vectors
for(n=0; n<nc; n++)
{

if(EdE==3 || EdE==6)
{
    if(a0_a[n]==0)
{
    cout << "\nmissing a0 for component " << n+1 << " in DATA BANK! \n";
    cout << "enter a0 value for component " << n+1 << ": ";
    cin >> a0[n];
}

    if(bCPA_a[n]==0)
{
    cout << "\nmissing bCPA for component " << n+1 << " in properties.csv! \n";
    cout << "enter bCPA value for component " << n+1 << ": ";
    cin >> bCPA[n];
}

    if(c1_a[n]==0)
{
    cout << "\nmissing c1 for component " << n+1 << " in properties.csv! \n";
    cout << "enter c1 value for component " << n+1 << ": ";
    cin >> c1[n];
}

    if(E_a[n]==0)
{
    cout << "\nmissing Epsilon(AB) for component " << n+1 << " in properties.csv! \n";
    cout << "enter Epsilon(AB) value for component " << n+1 << ": ";
    cin >> E[n];
}

    if(beta_a[n]==0)
{
    cout << "\nmissing beta(AB) for component " << n+1 << " in properties.csv! \n";
    cout << "enter beta(AB) value for component " << n+1 << ": ";
    cin >> beta[n];
}
}

Tc[n] = Tc_a[n];
Pc[n] = Pc_a[n];
omega[n] = omega_a[n];
r[n] = r_a[n];
q[n] = q_a[n];
q_prime[n] = q__a[n];
A[n] = A_a[n];
B[n] = B_a[n];
C[n] = C_a[n];
a0[n] = a0_a[n];
bCPA[n] = bCPA_a[n];
c1[n] = c1_a[n];
E[n] = E_a[n];
beta[n] = beta_a[n];
Tc_virtual[n] = Tc_v[n];
Pc_virtual[n] = Pc_v[n];
omega_virtual[n] = omega_v[n];
L_rg[n] = L_a[n];
phi_rg[n] = phi_a[n];
}

//--------------------------------------------------------------------------------
cout << "\n Choose the mixing rule: \n 1.Van der Waals \n 2.Van der Waals 2 (not working!) \n 3.Huron-Vidal" << endl;
cin >> MR;
output << "Mixing Rule = " << MR << endl;

if(MR!=1)
{
cout << "\n Choose the Excess Gibbs Energy Model: \n 1. UNIQUAC (not working!) \n 2.NRTL" << endl;
cin >> G_ex_model;
output << "Excess Gibbs Energy Model = " << G_ex_model << endl;

if(G_ex_model==2)
{
Aij(0,0) = 0;
Aij(1,1) = 0;
Aij(0,1) = 293.3380968; //Ethane+trifluoroethane 212.84K
Aij(1,0) = 286.330996;

alfa_NRTL(0,0) = 0;
alfa_NRTL(1,1) = 0;
alfa_NRTL(0,1) = 0.3;
alfa_NRTL(1,0) = alfa_NRTL(0,1);

cout << "\nA12= ";
cin >> Aij(0,1);
cout << "\nA12= ";
cin >> Aij(1,0);
cout << "\nalfa= ";
cin >> alfa_NRTL(0,1);
alfa_NRTL(1,0) = alfa_NRTL(0,1);
cout << "\nDij= ";
cin >> Dij; //Extra parameter when using CPA/HV-NRTL
cout << endl;
}

}

if(MR==1)
{
Dij = 1;
}

//Ideal gas constant
R = 0.08314462; // L.bar/K/mol
R = 0.000008314462; //MPa.m³/mol/K
//R = 8.3144598;

//Tolerances
tolZv = 0.0000001; //Erro para convergência de Zv
tolZl = 0.0000001; //Erro para convergência de Zl
tolSUMKx = 0.00001; //Erro para convergência no somatório de Kx
tolKx = 0.0001; //Erro para convergência de Kx
tolX = 0.0000001; //Fraction of non-associating sites tolerance
tolV = 0.000001; //Volume tolerance

Tr = T*Tc.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas
//Cálculo dos alfas
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);

//Updating EdE_parameters into vector
EdE_parameters = EdE_parameters_function(EdE);

//Calculating ai and bi
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);


cout << "\n The process is: \n 1.Isothermic \n 2.Isobaric" << endl;
cin >> process;

switch(process)
{
case 1: //Isothermic
cout << "\nDefine Temperature in K: ";
cin >> T;
output << "Defined Temperature = " << T << " K" << endl;
break;

case 2: //Isobaric
cout << "\nDefine Pressure in bar: ";
cin >> P;
//P = P*1e5; //from bar to Pa
output << "Defined Pressure = " << P << " kPa" << endl;
break;
}

cout << "\nUse binary interaction parameter? \n 1.Yes \n 2.No" << endl;
cin >> binary_interaction;
if(binary_interaction==1)
{
cout << "\nEnter the value for k12: ";
cin >> k12;
output << "kij = " << k12 << endl;
}

if(binary_interaction==2)
{
k12=0;
output << "kij = " << k12 << endl;
}

//SATURATION PRESSURE CALCULATION
switch(process)
{
case 1: //Isothermic
CT = C.array()+T;
//logPsat = A - (B*(CT.asDiagonal().inverse().diagonal().transpose())).diagonal();
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();
cout << "Psat bar= " << Psat << endl;
Psat = Psat.array()/10; //Converting from bar to MPa
cout << "Psat MPa= " << Psat << endl;
break;

case 2: //Isobaric
log10P = log10(P);
Alog10P = A.array()-log10P;
Tsat = (Alog10P.asDiagonal().inverse()*B)-C;
break;
}

//PRESSURE AND TEMPERATURE INITIAL GUESS
switch(process)
{
    case 1: //Isothermic
//Pressure initial guess
x(0) = 0.001;
x(1) = 1 - x(0);

if(mixture==1)
{
    x(0) = 0.999999;
    x(1) = 1 - x(0);
}

P = Psat.transpose()*x;
Pinit = P;
break;

    case 2: //Isobaric
//Temperature initial guess
x(0) = 0.001;
x(1) = 1 - x(0);

if(mixture==1)
{
    x(0) = 0.999999;
    x(1) = 1 - x(0);
}

Tb = Tsat.transpose()*x;
T = Tb;
Tinit = T;
CT = C.array()+T;
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();
/*
cout << "T = " << T << endl;
cout << "C = " << C << endl;
cout << "CT = " << CT << endl;
cout << "logPsat = " << logPsat << endl;
cout << "lnPsat = " << lnPsat << endl;
cout << "Isobaric Psat = " << Psat << endl;
*/
Pinit = P;
break;
}

if(EdE==3 || EdE==6)
{
    //Cross energy and volume of association calculation
    //The program asks the user for the combination rule in DELTA calculation
    cout << "\nChoose the combining rule:\n 1. CR-1 \n 2. CR-2 \n 3. CR-3 \n 4. CR-4 \n 5. ECR \n 6. CR-1 (+Solvation) (2B)" << endl;
    cin >> combining_rule;
    output << "Combining Rule = " << combining_rule << endl;

    if(combining_rule==6)
    {
        cout << "\nDefine BETCR: ";
        cin >> BETCR;
        cout << "\n";
        output << "BETCR = " << BETCR << endl;
    }

    if(combining_rule!=6)
    {
        BETCR=0.00;
    }

int nc4, i, j;
double E_component, beta_component;
int assoc_type[nc];
nc4 = 4*nc;



cout << "Consider the following choices for association schemes:\n 1. 1 \n 2. 2A \n 3. 2B \n 4. 3A \n 5. 3B \n 6. 4A \n 7. 4B \n 8. 4C" << endl;
for (j=0; j<nc; j++)
{
    cout <<   "Association scheme of component " << j+1 << ": ";
    cin >> assoc_type[j];
    output << "Association model component " << j+1 << " =" << assoc_type[j] << endl;
}

//Defining matrix for energy and volume of auto-association
for (j=0; j<nc; j++)
{

    E_component = E(j);
    beta_component = beta(j);
    E_i = energy_cross_association(assoc_type[j], E_component, R, T);
    beta_i = volume_cross_association(assoc_type[j], beta_component);

    E_n = energy_auto_association(assoc_type[j], E_component, R, T);
    beta_n = volume_auto_association(assoc_type[j], beta_component);

    E_col.block(4*j, 0, 4, 2) << E_i;
    beta_col.block(4*j, 0, 4, 2) = beta_i;

    E_auto.block(4*j,4*j,4,4) << E_n;
    beta_auto.block(4*j,4*j,4,4) << beta_n;

    E_row = E_col;
    beta_row = beta_col;
/*
    cout << "\nE_n = \n" << E_n << endl;
    cout << "\nE_auto = \n" << E_auto << endl;
    cout << "\nbeta_n = \n" << beta_n << endl;
    cout << "\nbeta_auto = \n" << beta_auto << endl;
*/
}
/*
for (i=0; i<nc; i++)
{
    E_component = E(i);
    beta_component = beta(i);
    E_i = energy_auto_association(assoc_type[i], E_component, R, T);
    beta_i = volume_auto_association(assoc_type[i], beta_component);

    for (j=0; j<nc; j++)
    {
        E_row.block(4*i,4*j,4,4) = E_i;
        beta_row.block(4*i,4*j,4,4) = beta_i;
    }
}
*/
}


//y initial guess
y = ((Psat*x.transpose()).diagonal()).array()/P;
yinit = y;


int iter_choice;
cout << "\nThe program calculates an automatic iteration for x1 going from 0.001 to 0.999, x2 = 1-x1" << endl;
cout << "Calculate a single point instead? \n 1.Yes \n 2.No" << endl;
cin >> iter_choice;

counter = 0;

output     << "-------------------------------------------------------------------------------------------------------------" << endl << endl;
output     << "Cálculos \n ----------------------------------------------------------------------------------------------------------" << endl;
output     << "x1 " << ";" << "y1 " << ";" << "T " << ";" << "P(kPa)" << ";" << "Vl" << ";"
           << "Vv" << ";" << "sumKx" << ";" << "counter" << ";" << "u_liquid1" << ";" << "u_vapor1" << ";"
           << "X1l" << ";" << "X1v" << ";" << "Zl" << ";" << "Zv" << ";" << "phi_liquid_1"
           << ";" << "phi_liquid_2" << ";" << "phi_vapor_1" << ";" << "phi_vapor_2" << ";"
           << "Vl_obj" << ";" << "Vv_obj" << ";" << "dP/dVl" << ";" << "dP/dVv" << ";"
           << "G_excess" << ";" << "P(bar)" << ";" << "density_l" << ";" << "density_v" << endl;

if(Renormalization==1)
{
    cout << "Renormalization density steps: ";
    cin >> n;

    cout << "everyone ok" << endl;

    if(EdE==1)
    {
        if(Tc_virtual[0] != 0)
        {
        Tr = T*Tc_virtual.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega_virtual, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, bCPA, EdE);
        }

        else
        {
        Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
        }

    }

    else
    {
    Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    ofstream Renorm;
    ofstream Not_splined;
    ofstream Envelope;
    ofstream dfnout;
    ofstream xpp_out;
    ofstream xp1_out;
    ofstream xp2_out;
    ofstream xpu_out;
    ofstream xpf_out;
    ofstream xpg_out;
    ofstream xpf1_out;
    ofstream xpu1_out;
    ofstream xpu2_out;
    //ofstream before_renorm;
    //ofstream after_renorm;

    Renorm.open("../Planilhas de análise/Renormalization.csv");
    Not_splined.open("../Planilhas de análise/Renormalization_not_splined.csv");
    dfnout.open("dfn_msa_out.csv");
    xpp_out.open("../Planilhas de análise/xpp.csv");
    xp1_out.open("../Planilhas de análise/xp1.csv");
    xp2_out.open("../Planilhas de análise/xp2.csv");
    xpu_out.open("../Planilhas de análise/xpu.csv");
    xpf_out.open("../Planilhas de análise/xpf.csv");
    xpg_out.open("../Planilhas de análise/xpg.csv");
    xpf1_out.open("../Planilhas de análise/xpf1.csv");
    xpu1_out.open("../Planilhas de análise/xpu1.csv");
    xpu2_out.open("../Planilhas de análise/xpu2.csv");
    //before_renorm.open("../Planilhas de análise/before_renorm.csv");
    //after_renorm.open("../Planilhas de análise/after_renorm.csv");
    Envelope.open("../Planilhas de análise/Envelope.csv");

    //Output files' headers
    Not_splined << "Density" << ";" << "f" << ";" << "f0" << ";" << "u" << ";" << "P" << ";" << "T" << endl;
    //before_renorm << "Density" << ";" << "f" << ";" << "flv+" << ";" << "flv-" << ";" << "f0a" << ";" << "T" << endl;
    Renorm << "Density" << ";" << "f" << ";" << "f0" << ";" << "u" << ";" << "P" << ";" << "u_0" << ";" << "P_0" << ";" << "T" << endl;
    Envelope << "rho" << ";" << "rho_l" << ";" << "P_l" << ";" << "rho_v" << ";" << "P_v" << endl;
    dfnout << "i" << ";" << "rho" << ";" << "f" << ";" << "delta_f" << endl;
    //after_renorm << "Density" << ";" << "f_before" << ";" << "f_after" << ";" << ";" << ";" << "T" << endl;

    int i, j, k, w, t;
    long double kB, L, L3, fi, K, rho, rho_plus, rho_minus, fl_plus, fl, fl_minus, fs_plus, fs, fs_minus;
    long double Gl, Gs, OMEGA, delta_f, f, f0, fl_old_plus, fl_old_minus, fs_old_plus, fs_old_minus, fl_old, fs_old;
    long double OMEGAs, OMEGAl, f_old, alfa_r, am, rho_max, bm, tolZ, rho2, var, f0_plus, f0_minus;
    long double width, suml, sums, m, fl_plus_old, fl_minus_old, fs_plus_old, fs_minus_old, f_original;
    long double Gl0, Gs0, Gln, Gsn, eGl0, eGs0, eGln, eGsn, phi_r, P_test_old, P_average_0;
    long double pmax_cond, P_max, P_min, P_average, P_test, test, P_l, u_l, P_v, u_v, pmin_cond, rho_v;
    double pi, eps, lambda, sigma, zeta_squared, NA;
    double sig1, sig2, sig12, l12, cnst;

    //MatrixXd Area(1000,1000);

    double Q_func, Kn, Ins, Inl, aminl, amins, al, as;
    double Area;
    int flag = 0;

    VectorXd L3_v(nc), L_v(nc), fi_v(nc), lnXX2(4*nc), f_assoc(nc), one_4c(4*nc, 2);

    VectorXd rhov(500), x_(500), fv_(500), X_plus(4*nc), X_minus(4*nc);


    std::vector<double> rho_vec_out(n), dP2dV2(n), dP_dV(n), P_vec(n), du_dV(n);
    std::vector<double> u_vec(n), u_vec_0(n), P_vec_0(n), f_vec_out(n), f0_vec_out(n);
    std::vector<double> rho1_vec_out(n), rho2_vec_out(n), u_vec_res1(n), u_vec_res2(n);
    //n = 5000;

    std::vector<double> rho_vec(n), f_vec(n), f_vec_res(n), u_vec1(n), f0_vec(n), P_vec1(n), Glv2(n), Gsv2(n);
    std::vector<double> rho1_vector(n), rho2_vector(n);
    std::vector<double> flvv(n), fsvv(n);
    VectorXd fl_old_p(n), fl_oldv(n), fl_old_m(n), fs_old_p(n), fs_oldv(n), fs_old_m(n), rho_vector2(n), f_after(n), f_before(n);
    VectorXd flv(n), fsv(n), fv(n), rho_vector(n), delta_fv(n), f_originalv(n), Glv(n), Gsv(n), argl(n), args(n);
    VectorXd P_vec_e(n), u_vec_e(n), rho_vec_out_e(n), f_vec_e(n), rho1_vec(n), rho2_vec(n);
    std::vector<double> V_vec(n), A_vec(n), f_env(n), rho_env(n), P_env(n), u_env(n);

    one_4c << 1, 0,
              1, 0,
              1, 0,
              1, 0,
              0, 1,
              0, 1,
              0, 1,
              0, 1;

    x << 0.999999, 0.000001;
    y << 0.999999, 0.000001;
    L_v << 7.40e-10, 7.50e-10;
    fi_v << 8.39, 8.98;

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    alfa_r = am;

    NA = 6.023e23;
    kB = 1.38064852e-25; // Boltzmann constant L.bar/K
    //kB = 1.38064852e-23; // J/K
    kB = R/NA;
    //kB = 1.38064852e-23; //Boltzmann constant J/K

    //L3_v = L_v.array().pow(3);
    //L3 = x.transpose()*L3_v;
    //L = pow(L3,(1/3)); //Cut-off Length
/*
    if(EdE==1)
    {
        cout << "L (dm) = ";
        cin >> L;
        cout << endl;

        cout << "phi = ";
        cin >> phi_r;
        cout << endl;
    }


    if(EdE==3)
    {
    if(cp[0]==47)
    {
        L = 5.51e-09;
        phi_r = 8.39;
    }//Initial shortest wavelength

    else
    {
        cout << "L = " << endl;
        cin >> L;
        cout << "phi_r = " << endl;
        cin >> phi_r;
    }
    }
    */

/*
if(EdE != 4)
{
        cout << "\nL (m) = ";
        cin >> L;
        cout << endl;

        cout << "\nphi = ";
        cin >> phi_r;
        cout << endl;
}
*/

    //fi = x.transpose()*fi_v; //Second crossover parameter phi

    rho_max = 0.99999/bm;

    //DIMENSIONLESS!!!************************************************************
    rho_max = 0.99999;

    t = 0;
    k = 0;
    w = 0;

    cout << "\nDefine final Temperature: ";
    cin >> final_T;

    cout << "\nDefine step size: ";
    cin >> step;

    int env_type;
    cout << "\nEnvelope type: \n1. Maxwell \n2. Area (not working) \n3.Newton \n4.Maxwell + Newton \n5.Faster Newton \n6.Newton seed\n";
    cin >> env_type;

    int critical_find;
    cout << "Type of envelope: \n1.Manual \n2.Automatic \n3.Find Tc fast \n4.Binary Mixture \n";
    cin >> critical_find;

    /*
    int estimation;
    cout << "Adjust to experimental data estimating L and phi?: \n1.Yes \n2.No \n";
    cin >> estimation;
    */

    int max_renorm;
    if(EdE==4) max_renorm = 6;
    if(EdE!=4) max_renorm = 9;

    int exponents;
    cout << "Calculate beta and delta critical exponents?: \n1.Yes \n2.No \n";
    cin >> exponents;


    init_T = T;
    Told = T;
    T = init_T;

    double T_original, step_original, final_T_original;
    T_original = T;
    step_original = step;
    final_T_original = final_T;
    int g;
    g = 0;

//while(T<(Tc[0]+10))

switch(critical_find)
{
    case 1: //Use designated final temperature and steps
while(T<=final_T)
{
    int p=0;

    if(step<0 && T<final_T/2) final_T = 0;

    //for(x(0)=0.000; x(0)<=1.000; x(0) = x(0) + 0.005)
    for(x(0)=0.000; x(0)<=1.000; x(0) = x(0) + 0.005)
    {
    if(x(0)==0) x(0)=0.000001;
    if(cp[0]==cp[1]) x(0) = 0.999999; //If both components are equal, then consider pure component
    cout << "x0 = " << x(0) << endl;

    x(1) = 1 - x(0);
    //RG parameters
    L = x(0)*pow(L_rg(0),3)+x(1)*pow(L_rg(1),3);
    L = pow(L,(1./3.));
    phi_r = x(0)*phi_rg(0) + x(1)*phi_rg(1);
    cout << "BEGIN==============================================================\n";
    //Recalculating Temperature-dependent parameters

    //DIMENSIONLESS!!!************************************************************
    T = T_original + g*step_original;
    final_T = final_T_original;
    step = step_original;

    cout << "T = " << T << endl;
    if(EdE==1)
    {
        if(Tc_virtual[0] != 0)
        {
        Tr = T*Tc_virtual.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega_virtual, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, bCPA, EdE);
        }

        else
        {
        Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
        }
    }

    else
    {
    Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    cout << "n = " << n << endl;

    if(EdE==4)
    {
        eps = 8.75;
        lambda = 1.65;
        sigma = 4.86e-10;
        L = 0.8e-4; //omega
        pi = 3.14159265359796;
        NA = 6.023e23;
        cnst = pi/6;
        l12 = 0;

        a = a_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0, T);
        b = b_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE, T);

        sig1 = pow((b(0)/cnst),(1.0/3.0));
        sig2 = pow((b(1)/cnst),(1.0/3.0));
        sig12 = 0.5*(sig1+sig2)*(1.0-l12);
        bm = (pi/6)*pow(sig12,3.0);
        am = pow((a(0)*a(1)),0.5)*(1.0-k12);
        /*
        cout << "sig1 = " << sig1 << endl;
        cout << "sig2 = " << sig2 << endl;
        cout << "sig12 = " << sig12 << endl;
        cout << "b(0) = " << b(0) << endl;
        cout << "b(1) = " << b(1) << endl;
        cout << "cnst = " << cnst << endl;
        cout << "bm = " << bm << endl;
        cout << "am = " << am << endl;
        */

        zeta_squared = pow(lambda,2)*pow(sigma,2)/5;

        rho_max = 6/(pi*pow(sigma,3))/NA;

        //DIMENSIONLESS
        rho_max = rho_max*bm;

        phi_r = zeta_squared*10;
    }

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    w = 0;

    //DIMENSIONLESS!!!************************************************************
    T = T*bm*R/am;
    final_T = final_T_original*bm*R/am;
    step = step_original*bm*R/am;

    if(EdE==4)
    {
        T = T/bm;
        final_T = final_T/bm;
        step = step/bm;
    }


    //====================================================================
    cout << "Before renormalization =======================================\n";

if(r_type==1)
{
    //Calcular vetor de f em f0 com um cálculo
    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm; //adimensionalization

    for(k=0; k<n; k++)
    {
    rho_vec[k] = double(k)/n/bm;
    rho_vector(k) = double(k)/n/bm;
    rho_vec[0] = 1e-5;
    rho_vector(0) = 1e-5;

    //NON BONDED FRACTION DIMENSIONAL**************************
    T = T/bm/R*am;

    if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_vector(k), deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

    T = T*bm*R/am; //adimensionalization
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    rho_vec[k] = rho_vec[k]*bm; //adimensionalization
    rho_vector(k) = rho_vector(k)*bm; //adimensionalization
    rho_vec[0] = rho_vec[0]*bm; //adimensionalization
    rho_vector(0) = rho_vector(0)*bm; //adimensionalization

    fv(k) = helmholtz_repulsive(EdE, R, T, rho_vector(k), am, bm, X, x, sigma, eps, kB);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << endl;

    //DIMENSIONLESS!!!************************************************************
    ////////if(EdE != 4) fv(k) = fv(k) + 0.5*am*rho_vector(k)*rho_vector(k);

    if(EdE != 4) fv(k) = fv(k) + 0.5*rho_vector(k)*rho_vector(k);
    //DIMENSIONLESS!!!************************************************************

    f_originalv(k) = fv(k);

    //DIMENSIONLESS!!!************************************************************
    /////////////if(EdE != 4) f_originalv(k) = fv(k) - 0.5*am*rho_vector(k)*rho_vector(k);


    if(EdE != 4) f_originalv(k) = fv(k) - 0.5*rho_vector(k)*rho_vector(k);
    //DIMENSIONLESS!!!************************************************************


    f0_vec[k] = f_originalv(k);
    f_before(k) = fv(k);
    rho = rho + rho_max/n;
    }

    //Iteração principal - o n
    for(i=1; i<9; i++) //i começando em 1
    {
        //Calcular K
        Kn = kB*T/((pow(2,3*i))*pow(L,3));

        //DIMENSIONLESS!!!************************************************************
        Kn = T/pow(2,3*i)/(pow(L,3)/bm*6.023e23);


        //DIMENSIONLESS FOR MSA!!!!!
        if(EdE==4)
        {
        T = T/R*am;

        Kn = R*T/((pow(2,3*i))*L);
        Kn = Kn/am*bm;

        T = T*R/am;
        }

        //Kn = kB*T*NA/((pow(2,3*i))*L*b(0));

        //Preparar vetores f_l e f_s
        rho = 1e-6;

        if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

        //DIMENSIONLESS!!!************************************************************
        rho = 1e-6*bm;

        for(w=0; w<n; w++)
        {
        flv(w) = helmholtz_recursion_long(EdE, fv(w), rho_vector(w), am, bm, R, T);
        fsv(w) = helmholtz_recursion_short(EdE, fv(w), rho_vector(w), am, bm, i, L, phi_r, sr_type, R, T);
        rho = rho + rho_max/n;
        }

        //cout << "flv before 201 / 201-100 = " << flv(201) << " / " << flv(201-100) << endl;
        //cout << "fsv before 201 / 201-100 = " << fsv(201) << " / " << fsv(201-100) << endl;

            //Iteração 2 - calcular os valores para f no i atual
            rho = 1e-6;

            //DIMENSIONLESS!!!************************************************************
            rho = 1e-6*bm;

            w = 0;
            width = rho_max/n;

            for(w=0; w<n; w++)
            {
                if(EdE != 4) delta_fv(w) = df_calculation(w,n,Iteration,width,Kn,rho_vec,flv,fsv);
                if(EdE == 4)
                {
                    if(w<n/2) delta_fv(w) = df_calculation(w,n,Iteration,width,Kn,rho_vec,flv,fsv);
                    else delta_fv(w) = 1e-15;
                }
            }


        //Calcular o novo vetor de f, ajustando com o vetor de delta_f
        fv.array() = fv.array() + delta_fv.array();
        cout << "i = " << i << " / Kn = " << Kn*am/bm << " / f 60 = " << fv(60)*am/bm << " / delta 60 --> "
             << delta_fv(60)*am/bm << endl;

        //cin >> stop;
        /*
        for(w=0; w<n; w++)
        {
            dfnout << i << ";" << rho_vector(w)/bm << ";" << fv(w)*am/bm << ";" << delta_fv(w)*am/bm << endl;
        }
        */
    }

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    for(w=0; w<n; w++)
    {
    //DIMENSIONLESS!!!************************************************************
    /////////////if(EdE != 4) fv(w) = fv(w) - 0.5*am*rho_vector(w)*rho_vector(w);


    if(EdE != 4) fv(w) = fv(w) - 0.5*rho_vector(w)*rho_vector(w);
    //DIMENSIONLESS!!!************************************************************

    f_vec[w] = fv(w);

    //DIMENSIONLESS!!!************************************************************
    if(EdE != 4)
    {
    f_vec[w] = fv(w)*am/bm/bm; //Dimensionalization
    f0_vec[w] = f0_vec[w]*am/bm/bm;
    rho_vec[w] = rho_vec[w]/bm; //Dimensionalization
    }

    if(EdE == 4)
    {
    f_vec[w] = fv(w)*am/bm;
    f0_vec[w] = f0_vec[w]*am/bm;
    rho_vec[w] = rho_vec[w]/bm;
    }


    //cout << "rho = " << rho_vector(w) << "  //  f = " << fv(w) << endl;
    Not_splined << std::fixed << std::setprecision(15) << rho_vector(w) << ";" << fv(w) << ";" << f_originalv(w) << ";"
                << T << endl;
    rho = rho + rho_max/n;
    }
}
    //====================================================================

rho = 1e-6;
w=0;
double rho1v[n], rho2v[n];

for(w=0; w<n; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/n/bm;
    rho_vec_out[0] = 1e-6;

    rho1_vec[w] = double(w)/n/b(0);
    rho2_vec[w] = double(w)/n/b(1);
    rho1_vec_out[w] = double(w)/n/b(0);
    rho2_vec_out[w] = double(w)/n/b(1);
    rho1v[w] = rho1_vec_out[w];
    rho2v[w] = rho2_vec_out[w];
    rho1_vector[w] = rho1_vec[w];
    rho2_vector[w] = rho2_vec[w];

    //rho1_vec[w] = rho_vec_out[w]*x(0);
    //rho2_vec[w] = rho_vec_out[w]*x(1);
    //rho1_vec_out[w] = rho_vec_out[w]*x(0);
    //rho2_vec_out[w] = rho_vec_out[w]*x(1);

    //cout << "densities m 1 2: " << rho_vec_out[w] << " " << rho1_vec_out[w] << " " << rho2_vec_out[w] << endl;

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

    if(EdE!=4) T = T*am/bm/R;
    if(EdE==4) T = T*am/R;

//Subtract ideal gas contribution before cubic spline
for(w=0; w<n; w++)
{
    f_vec[w] = f_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f0_vec[w] = f0_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f_vec_res[w] = f_vec[w];
    //f_mat[w][k2] = f_vec[w];
}

f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);

u_vec_res1 = cspline_deriv1_vec(rho1_vector, f_vec, rho1_vec_out);
u_vec_res2 = cspline_deriv1_vec(rho2_vector, f_vec, rho2_vec_out);
//Add ideal gas contribution after cubic spline
for(w=0; w<n; w++)
{
    f_vec[w] = f_vec[w] + rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f_vec_out[w] = f_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    f0_vec_out[w] = f0_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);

    //SHOULDN'T IT BE xi FOR MIXTURES INSTEAD OF RHO???
    u_vec[w] = u_vec[w] + R*T*(log(rho_vec_out[w]));
    u_vec_0[w] = u_vec_0[w] + R*T*(log(rho_vec_out[w]));
}


    if(EdE!=4) T = T/am*bm*R;
    if(EdE==4) T = T/am*R;

for(i=0; i<n; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b1;
    a = (P_vec[n/2+5]-P_vec[n/2-5])/(rho_vec_out[n/2+5]-rho_vec_out[n/2-5]);
    b1 = P_vec[n/2+5]-a*rho_vec_out[n/2+5];
    P_vec[n/2-4] = a*rho_vec_out[n/2-4] + b1;
    P_vec[n/2-3] = a*rho_vec_out[n/2-3] + b1;
    P_vec[n/2-2] = a*rho_vec_out[n/2-2] + b1;
    P_vec[n/2-1] = a*rho_vec_out[n/2-1] + b1;
    P_vec[n/2] = a*rho_vec_out[n/2] + b1;
    P_vec[n/2+1] = a*rho_vec_out[n/2+1] + b1;
    P_vec[n/2+2] = a*rho_vec_out[n/2+2] + b1;
    P_vec[n/2+3] = a*rho_vec_out[n/2+3] + b1;
    P_vec[n/2+4] = a*rho_vec_out[n/2+4] + b1;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<n; i++)
{
 //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
 //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;

        //DIMENSIONLESS!!!************************************************************
 if(EdE!=4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/bm/R << endl;
 }

 if(EdE==4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i]<< ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/R << endl;
 }

}

//=============================================================

if(EdE==4) cout << "T = " << T*am/R << endl;
if(EdE!=4) cout << "T = " << T*am/bm/R << endl;
cout << "=======================================\n" << endl;

    //****************************START WRITING FILE WITH f function of mole fraction of density*****************
    if(p<1)
    {
        xpp_out << ";"; //Cell A1 clear
        xp1_out << ";"; //Cell A1 clear
        xp2_out << ";"; //Cell A1 clear
        xpu_out << ";"; //Cell A1 clear
        xpf_out << ";"; //Cell A1 clear
        xpf1_out << ";"; //Cell A1 clear
        xpg_out << ";"; //Cell A1 clear
        xpu1_out << ";"; //Cell A1 clear
        xpu2_out << ";"; //Cell A1 clear

        //xpu1_out << 0 << ";"; //Cell A2 0
        //xpu2_out << 0 << ";"; //Cell A2 0

        //Write density values on first line
        for(int i=0; i<n; i++)
        {
        xpp_out << rho_vec_out[i]*bm << ";";//Handle P with mole fraction and density
        xp1_out << rho1_vec_out[i]*b(0) << ";";//Handle P with mole fraction and density
        xp2_out << rho2_vec_out[i]*b(1) << ";";//Handle P with mole fraction and density
        xpu_out << rho_vec_out[i]*bm << ";";//Handle u with mole fraction and density
        xpf_out << rho_vec_out[i]*bm << ";";//Handle u with mole fraction and density
        xpg_out << rho_vec_out[i]*bm << ";";//Handle u with mole fraction and density
        //xpf1_out << (1/(rho_vec_out[n-1-i]*bm)) << ";";//Handle u with mole fraction and density
        xpf1_out << (1/(rho_vec_out[n-1-i]*bm)) << ";";//Handle u with mole fraction and density
        xpu1_out << rho1_vec_out[i]*b(0) << ";";//Handle u with mole fraction and density
        xpu2_out << rho2_vec_out[i]*b(1) << ";";//Handle u with mole fraction and density
        }

        xpp_out << endl; //Jump to next line, start mole fractions P and u
        xp1_out << endl; //Jump to next line, start mole fractions P and u
        xp2_out << endl; //Jump to next line, start mole fractions P and u
        xpu_out << endl; //Jump to next line, start mole fractions P and u
        xpf_out << endl; //Jump to next line, start mole fractions P and u
        xpg_out << endl; //Jump to next line, start mole fractions P and u
        xpf1_out << endl; //Jump to next line, start mole fractions P and u
        xpu1_out << endl; //Jump to next line, start mole fractions P and u
        xpu2_out << endl; //Jump to next line, start mole fractions P and u
    }


    xpp_out << x(0) << ";"; //Write mole fraction on first column
    xp1_out << x(0) << ";"; //Write mole fraction on first column
    xp2_out << x(0) << ";"; //Write mole fraction on first column
    xpu_out << x(0) << ";"; //Write mole fraction on first column
    xpf_out << x(0) << ";"; //Write mole fraction on first column
    xpg_out << x(0) << ";"; //Write mole fraction on first column
    xpf1_out << x(0) << ";"; //Write mole fraction on first column
    xpu1_out << x(0) << ";"; //Write mole fraction on first column
    xpu2_out << x(0) << ";"; //Write mole fraction on first column DOUBT

    for(i=0; i<n; i++)
    {
        xpp_out << P_vec[i] << ";";//Handle P with mole fraction and density
        xp1_out << P_vec[i] << ";";//Handle P with mole fraction and density
        xp2_out << P_vec[i] << ";";//Handle P with mole fraction and density
        xpu_out << u_vec[i] << ";";//Handle u with mole fraction and density
        xpf_out << f_vec_res[i] << ";";//Handle u with mole fraction and density
        //xpf1_out << f_vec[n-1-i]*bm*bm/am << ";";//dimensionless
        xpf1_out << f_vec[i] << ";";
        xpg_out << (f_vec[n-1-i] + P_vec[n-1-i])/(rho_vec_out[n-1-i]) << ";";
        xpu1_out << u_vec_res1[i] << ";";//Handle u with mole fraction and density
        xpu2_out << u_vec_res2[i] << ";";//Handle u with mole fraction and density
    }

        xpp_out << bm << endl; //Jump to next line, next mole fraction
        xp1_out << b(0) << endl; //Jump to next line, next mole fraction
        xp2_out << b(1) << endl; //Jump to next line, next mole fraction
        xpu_out << bm << endl; //Jump to next line, next mole fraction
        xpf_out << bm << endl; //Jump to next line, next mole fraction
        xpg_out << bm << endl; //Jump to next line, next mole fraction
        xpf1_out << bm << endl; //Jump to next line, next mole fraction
        xpu1_out << b(0) << endl; //Jump to next line, next mole fraction
        xpu2_out << b(1) << endl; //Jump to next line, next mole fraction
    //****************************END WRITING FILE WITH f function of mole fraction of density*****************

k++;
p++;

    if(x(0)==0.000001) x(0)=x(0)-0.000001;
    } //end for x loop
    T = T + step;
    g++;
} //end while T loop

    if(cp[0] != cp[1])
    {
/*
    //Begin writing f with rho1 and rho2 to find critical point**************************************
    xpf1_out << ";"; //Cell A1 clear
    for(int i=0; i<n; i++)
    {
    xpf1_out << rho2_vec_out[i] << ";";//Write density2 on first line
    }
    xpf1_out << endl; //Jump to next line, start mole fractions P and u

    for(j=0; j<n; j++)
    {
        xpf1_out << rho1_vec_out[j] << ";"; //Write density1 on first column
        for(i=0; i<n; i++)
        {
        xpf1_out << f_vec[i] << ";";//Handle u with mole fraction and density
        }
        xpf1_out << bm << endl; //End line with bm and jump to next line, next density1
    }
    //End writing f with rho1 and rho2 to find critical point**************************************
*/
        d2P = new double *[200];
        for(int k = 0; k <200; k++)
            d2P[k] = new double[n];

        d2P1 = new double *[200];
        for(int k = 0; k <200; k++)
            d2P1[k] = new double[n];

        d2P2 = new double *[200];
        for(int k = 0; k <200; k++)
            d2P2[k] = new double[n];

        d2u1 = new double *[200];
        for(int k = 0; k <200; k++)
            d2u1[k] = new double[n];

        d2u2 = new double *[200];
        for(int k = 0; k <200; k++)
            d2u2[k] = new double[n];

        Pmat = new double *[200];
        for(int k = 0; k <200; k++)
            Pmat[k] = new double[n];

        P1mat = new double *[200];
        for(int k = 0; k <200; k++)
            P1mat[k] = new double[n];

        P2mat = new double *[200];
        for(int k = 0; k <200; k++)
            P2mat[k] = new double[n];

        umat = new double *[200];
        for(int k = 0; k <200; k++)
            umat[k] = new double[n];

        u1mat = new double *[200];
        for(int k = 0; k <200; k++)
            u1mat[k] = new double[n];

        u2mat = new double *[200];
        for(int k = 0; k <200; k++)
            u2mat[k] = new double[n];

        Wmat = new double *[200];
        for(int k = 0; k <200; k++)
            Wmat[k] = new double[n];

        Dmat = new double *[200];
        for(int k = 0; k <200; k++)
            Dmat[k] = new double[n];

        rho_rv.resize(n);

        cout << "before hicks young" << endl;

        //hicks_young_critical(Wmat,Dmat,n);

        cout << "after hicks young" << endl;

        d2Pgen(d2P,n);
        d2P1gen(d2P1,n);
        d2P2gen(d2P2,n);
        cout << "d2P ok" << endl;
        d2u1gen(d2u1,n);
        d2u2gen(d2u2,n);
        //new_d2u1gen(d2u1,b(0),b(1),n);
        //new_d2u2gen(d2u2,b(0),b(1),n);
        cout << "d2u1 and d2u2 ok" << endl;
        renorm_mat_reader(Pmat,umat,n);
        cout << "Pmat umat ok" << endl;
        renorm_mat_reader(P1mat,P2mat,n);
        cout << "P1mat P2mat ok" << endl;
        //renorm_uu_reader(u1mat,u2mat,n);
        //cout << "u1mat u2mat ok" << endl;
        rho_rv = renorm_rhovec(n);
        cout << "rho_rv ok" << endl;
        x_rv = renorm_xvec(n);
        cout << "x_rv ok" << endl;
        d2ugen(d2u1,d2u2,u1mat,u2mat,n);
        cout << "new u ok" << endl;
        T = T_original;
        cout << "T = " << T << endl;
        Renormalization=2;
        if(EdE==1) EdE = 5;
        if(EdE==3) EdE = 6;
    }
    else envelope_tracer(1e-5,env_type,n);

    break;

 //===============================================================================================================//
 //
 //                                       AUTOMATIC CRITICAL POINT BELOW
 //
 //==============================================================================================================//

    case 2: //Binary mixture using MSA

    double delta_fm;

    double **fm;
    double **f_mat_res;
    double **Pmat;
    double **f_originalm;
    double **f0_mat;
    double **rho_mat;
    double **flm;
    double **fsm;
    double **umat;
    double **fr_mat;
    double **fmm;
    double **fn1;
    double **fn2;
    double **fn11;
    double **fn22;
    double **P_mat;

    fm = new double *[n];
    f_mat_res = new double *[n];
    Pmat = new double *[n];
    f_originalm = new double *[n];
    f0_mat = new double *[n];
    rho_mat = new double *[n];
    flm = new double *[n];
    fsm = new double *[n];
    umat = new double *[n];
    fr_mat = new double *[n];
    fmm = new double *[n];
    fn1 = new double *[n];
    fn2 = new double *[n];
    fn11 = new double *[n];
    fn22 = new double *[n];
    P_mat = new double *[n];

    for(int k=0; k<n; k++)
    {
    fm[k] = new double[n];
    f_mat_res[k] = new double[n];
    Pmat[k] = new double[n];
    f_originalm[k] = new double[n];
    f0_mat[k] = new double[n];
    rho_mat[k] = new double[n];
    flm[k] = new double[n];
    fsm[k] = new double[n];
    umat[k] = new double[n];
    fr_mat[k] = new double[n];
    fmm[k] = new double[n];
    fn1[k] = new double[n];
    fn2[k] = new double[n];
    fn11[k] = new double[n];
    fn22[k] = new double[n];
    P_mat[k] = new double[n];
    }


    while(T<=final_T)
{
    int p=0;
    int k1=0;

    //for(x(0)=0.000; x(0)<=1.000; x(0) = x(0) + 0.005)
    //{
    //if(x(0)==0) x(0)=0.000001;
    //cout << "x0 = " << x(0) << endl;

    //if(cp[0]==cp[1]) x(0) = 1; //If both components are equal, then consider pure component

    //x(1) = 1 - x(0);
    //RG parameters
    //L = x(0)*pow(L_rg(0),3)+x(1)*pow(L_rg(1),3);
    //L = pow(L,(1./3.));
    //phi_r = x(0)*phi_rg(0) + x(1)*phi_rg(1);
    cout << "BEGIN==============================================================\n";
    //Recalculating Temperature-dependent parameters

    //DIMENSIONLESS!!!************************************************************
    T = T_original + g*step_original;
    final_T = final_T_original;
    step = step_original;

    cout << "T = " << T << endl;
    if(EdE==1)
    {
        if(Tc_virtual[0] != 0)
        {
        Tr = T*Tc_virtual.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega_virtual, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, bCPA, EdE);
        }

        else
        {
        Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
        }
    }

    else
    {
    Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    cout << "n = " << n << endl;

    if(EdE==4)
    {
        eps = 8.75;
        lambda = 1.65;
        sigma = 4.86e-10;
        L = 0.8e-4; //omega
        pi = 3.14159265359796;
        NA = 6.023e23;
        cnst = pi/6;
        l12 = 0;

        a = a_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0, T);
        b = b_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE, T);

        sig1 = pow((b(0)/cnst),(1.0/3.0));
        sig2 = pow((b(1)/cnst),(1.0/3.0));
        sig12 = 0.5*(sig1+sig2)*(1.0-l12);
        bm = (pi/6)*pow(sig12,3.0);
        am = pow((a(0)*a(1)),0.5)*(1.0-k12);
        /*
        cout << "sig1 = " << sig1 << endl;
        cout << "sig2 = " << sig2 << endl;
        cout << "sig12 = " << sig12 << endl;
        cout << "b(0) = " << b(0) << endl;
        cout << "b(1) = " << b(1) << endl;
        cout << "cnst = " << cnst << endl;
        cout << "bm = " << bm << endl;
        cout << "am = " << am << endl;
        */


        zeta_squared = pow(lambda,2)*pow(sigma,2)/5;

        rho_max = 6/(pi*pow(sigma,3))/NA;

        //DIMENSIONLESS
        rho_max = rho_max*bm;

        phi_r = zeta_squared*10;
    }

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    w = 0;

    //DIMENSIONLESS!!!************************************************************
    //T = T*bm*R/am;
    //final_T = final_T_original*bm*R/am;
    //step = step_original*bm*R/am;

    if(EdE==4)
    {
        T = T/bm;
        final_T = final_T/bm;
        step = step/bm;
    }


    //====================================================================
    cout << "Before renormalization =======================================\n";

if(r_type==1)
{
    //Calcular vetor de f em f0 com um cálculo
    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm; //adimensionalization

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
    //MSA
    rho1_vector[i] = double(i)/n/b(0);
    rho2_vector[j] = double(j)/n/b(1);
    rho1_vec(i) = double(i)/n/b(0);
    rho2_vec(j) = double(j)/n/b(1);
    rho_mat[i][j] = x(0)*rho1_vector[i]+x(1)*rho2_vector[j];

    //NON BONDED FRACTION DIMENSIONAL**************************
    //T = T/bm/R*am;

    if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_mat[i][j], deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

    //T = T*bm*R/am; //adimensionalization
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    //rho_vec[k] = rho_vec[k]*bm; //adimensionalization
    //rho_vector(k) = rho_vector(k)*bm; //adimensionalization
    //rho_vec[0] = rho_vec[0]*bm; //adimensionalization
    //rho_vector(0) = rho_vector(0)*bm; //adimensionalization
    double ytot = rho1_vector[k]*b(0)+rho2_vector[k]*b(1);
    rho1_vec(k) = rho1_vec(k)*b(0);
    rho2_vec(k) = rho2_vec(k)*b(1);

    //T = T/bm/R*am; //dimensional
    fm[i][j] = helmholtz_repulsive_m(EdE, R, T, rho1_vector[i], rho2_vector[j], am, b(0), b(1), bm, X, x, sigma, eps, kB);
    if(ytot<0.9 && EdE==4) fm[i][j] = pow(1,30);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << endl;

    //DIMENSIONLESS!!!************************************************************
    //////////if(EdE != 4) fv(k) = fv(k) + 0.5*am*rho_vector(k)*rho_vector(k);

    if(EdE != 4) fm[i][j] = fm[i][j] + 0.5*rho_mat[i][j]*rho_mat[i][j];
    //DIMENSIONLESS!!!************************************************************

    f_originalm[i][j] = fm[i][j];

    //DIMENSIONLESS!!!************************************************************
    /////////////if(EdE != 4) f_originalv(k) = fv(k) - 0.5*am*rho_vector(k)*rho_vector(k);


    if(EdE != 4) f_originalm[i][j] = fm[i][j] - 0.5*rho_mat[i][j]*rho_mat[i][j];
    //DIMENSIONLESS!!!************************************************************


    f0_mat[i][j] = f_originalm[i][j];
    //f_before(k) = fv(k);
    //rho = rho + rho_max/n;
        }
    }

    //Iteração principal - o n
    for(i=1; i<9; i++) //i começando em 1
    {
        //Calcular K
        Kn = kB*T/((pow(2,3*i))*pow(L,3));

        //DIMENSIONLESS!!!************************************************************
        //Kn = T/pow(2,3*i)/(pow(L,3)/bm*6.023e23);


        //DIMENSIONLESS FOR MSA!!!!!
        if(EdE==4)
        {
        T = T/R*am;

        Kn = R*T/((pow(2,3*i))*L);
        Kn = Kn/am*bm;

        T = T*R/am;
        }

        //Kn = kB*T*NA/((pow(2,3*i))*L*b(0));

        //Preparar vetores f_l e f_s
        rho = 1e-6;

        if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

        //DIMENSIONLESS!!!************************************************************
        rho = 1e-6*bm;

        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
            flm[i][j] = helmholtz_recursion_long_m(EdE, fm[i][j], rho1_vector[i], rho2_vector[j], am, b(0), b(1), bm, R, T, rho_mat[i][j]);
            fsm[i][j] = helmholtz_recursion_short_m(EdE, fm[i][j], rho1_vector[i], rho2_vector[j], am, b(0), b(1), bm, i, L,
                                                    phi_r, sr_type, R, T, rho_mat[i][j]);
            }
        }

        //cout << "flv before 201 / 201-100 = " << flv(201) << " / " << flv(201-100) << endl;
        //cout << "fsv before 201 / 201-100 = " << fsv(201) << " / " << fsv(201-100) << endl;

            w = 0;
            width = rho_max/n;
            for(int i=0; i<n; i++)
            {
            for(int j=0; j<n; j++)
            {
                //if(EdE != 4) delta_fv(w) = df_calculation(w,n,Iteration,width,Kn,rho_vec,flv,fsv);
                if(EdE != 4) delta_fm = 1e-15;
                if(EdE == 4)
                {
                    //if(w<n/2) delta_fm(w) = dfm_calculation(i,j,n,Iteration,width,Kn,rho1_vec,rho2_vec,flm,fsm);
                    //else delta_fm(w) = 1e-15;
                    delta_fm = 1e-15;
                }
                fm[i][j] = fm[i][j]+delta_fm;
                if(i%50==0 && j%50==0) cout << "i,j= " << i << " " << j << " " << fm[i][j] << " " << delta_fm << endl;
            }
            }


        //Calcular o novo vetor de f, ajustando com o vetor de delta_f
        //fv.array() = fv.array() + delta_fv.array();
        //cout << "i = " << i << " / Kn = " << Kn*am/bm << " / f 60 = " << fv(60)*am/bm << " / delta 60 --> "
        //     << delta_fv(60)*am/bm << endl;
        //cin >> stop;
        /*
        for(w=0; w<n; w++)
        {
            dfnout << i << ";" << rho_vector(w)/bm << ";" << fv(w)*am/bm << ";" << delta_fv(w)*am/bm << endl;
        }
        */


    }

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
    //DIMENSIONLESS!!!************************************************************
    /////////////if(EdE != 4) fv(w) = fv(w) - 0.5*am*rho_vector(w)*rho_vector(w);


    if(EdE != 4) fm[i][j] = fm[i][j] - 0.5*rho_mat[i][j]*rho_mat[i][j];
    //DIMENSIONLESS!!!************************************************************

    //f_vec[w] = fv(w);

    //DIMENSIONLESS!!!************************************************************
    if(EdE != 4)
    {
    //f_mat[i][j] = f_mat(w)*am/bm/bm; //Dimensionalization
    //f0_vec[w] = f0_vec[w]*am/bm/bm;
    //rho_mat[i][j] = rho1_vector[i]/b(0)+; //Dimensionalization
    }

    if(EdE == 4)
    {
    //f_mat[w] = fm(w)*am/bm;

    //f0_vec[w] = f0_vec[w]*am/bm;
    //rho_vec[w] = rho_vec[w]/bm;
    //rho1_vec[w] = rho1_vec[w]/b(0);
    //rho2_vec[w] = rho2_vec[w]/b(1);
    }


    //cout << "rho = " << rho_vector(w) << "  //  f = " << fv(w) << endl;
    Not_splined << std::fixed << std::setprecision(15) << rho_vector(w) << ";" << fv(w) << ";" << f_originalv(w) << ";"
                << T << endl;
    rho = rho + rho_max/n;
        }
    }
}
    //====================================================================



rho = 1e-6;
w=0;

for(w=0; w<n; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/n/bm;
    rho_vec_out[0] = 1e-6;

    rho1_vec[w] = double(w)/n/b(0);
    rho2_vec[w] = double(w)/n/b(1);
    rho1_vec_out[w] = double(w)/n/b(0);
    rho2_vec_out[w] = double(w)/n/b(1);

    //rho1_vec[w] = rho_vec_out[w]*x(0);
    //rho2_vec[w] = rho_vec_out[w]*x(1);
    //rho1_vec_out[w] = rho_vec_out[w]*x(0);
    //rho2_vec_out[w] = rho_vec_out[w]*x(1);

    //cout << "densities m 1 2: " << rho_vec_out[w] << " " << rho1_vec_out[w] << " " << rho2_vec_out[w] << endl;

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

    //if(EdE!=4) T = T*am/bm/R;
    if(EdE==4) T = T*am/R;

//Subtract ideal gas contribution before cubic spline
for(w=0; w<n; w++)
{
    f_vec[w] = f_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f0_vec[w] = f0_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
}

//Spline fit free energy
double rho1v[n], rho2v[n];
    for(int i=0; i<n; i++)
    {
    rho1v[i]=rho1_vector[i];
    rho2v[i]=rho2_vector[i];
    }

for(int i=0; i<n; i++)
{
    for(int j=0; j<n; j++)
    {
    //f_mat[i][j] = fm[i][j];
    }
}

for(int i=0; i<n; i++)
{
    for(int j=0; j<n; j++)
    {
    fr_mat[i][j] = fm[i][j] - f0_mat[i][j];
    fmm[j][i] = fm[i][j];
    rho_mat[i][j] = x(0)*rho1_vec[i]+x(1)*rho2_vec[j];
    f_mat_res[i][j] = fm[i][j] - rho_mat[i][j]*R*T*(log(rho_mat[i][j])-1);
    }
}

splie2(rho1v,rho2v,n,n,fm,fn22);
splie2(rho2v,rho1v,n,n,fmm,fn11);
forward_row_deriv(rho1v,rho2v,n,n,f_mat_res,fn1);
forward_col_deriv(rho1v,rho2v,n,n,f_mat_res,fn2);

double tmp;

for(int i=0; i<n-1; i++)
{
    for(int j=i+1; j<n; j++)
    {
     tmp = fn11[i][j];
     fn11[i][j] = fn11[j][i];
     fn11[j][i] = tmp;
    }
}

for(int i=0; i<n; i++)
{
    for(int j=0; j<n; j++)
    {
     P_mat[i][j] = rho1_vector[i]*fn1[i][j] + rho2_vector[j]*fn2[i][j] - fm[i][j];
    }
}

//f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
//f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

//u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
//u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);

//u_vec_res1 = cspline_deriv1_vec(rho1_vec, f_vec, rho1_vec_out);
//u_vec_res2 = cspline_deriv1_vec(rho2_vec, f_vec, rho2_vec_out);

//u_vec_res1 = cspline_deriv1_vec(rho1_vector, f_vec, rho1_vec_out);
//u_vec_res2 = cspline_deriv1_vec(rho2_vector, f_vec, rho2_vec_out);

//Add ideal gas contribution after cubic spline
for(w=0; w<n; w++)
{
    //f_vec_out[w] = f_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    //f0_vec_out[w] = f0_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);

    //SHOULDN'T IT BE xi FOR MIXTURES INSTEAD OF RHO???
    //u_vec[w] = u_vec[w] + R*T*(log(rho_vec_out[w]));
    //u_vec_0[w] = u_vec_0[w] + R*T*(log(rho_vec_out[w]));
}
    //if(EdE!=4) T = T/am*bm*R;
    if(EdE==4) T = T/am*R;

for(i=0; i<n; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b1;
    a = (P_vec[n/2+5]-P_vec[n/2-5])/(rho_vec_out[n/2+5]-rho_vec_out[n/2-5]);
    b1 = P_vec[n/2+5]-a*rho_vec_out[n/2+5];
    P_vec[n/2-4] = a*rho_vec_out[n/2-4] + b1;
    P_vec[n/2-3] = a*rho_vec_out[n/2-3] + b1;
    P_vec[n/2-2] = a*rho_vec_out[n/2-2] + b1;
    P_vec[n/2-1] = a*rho_vec_out[n/2-1] + b1;
    P_vec[n/2] = a*rho_vec_out[n/2] + b1;
    P_vec[n/2+1] = a*rho_vec_out[n/2+1] + b1;
    P_vec[n/2+2] = a*rho_vec_out[n/2+2] + b1;
    P_vec[n/2+3] = a*rho_vec_out[n/2+3] + b1;
    P_vec[n/2+4] = a*rho_vec_out[n/2+4] + b1;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<n; i++)
{
// Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
//        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;

        //DIMENSIONLESS!!!************************************************************
 if(EdE!=4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/bm/R << endl;
 }

 if(EdE==4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i]<< ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/R << endl;
 }

}

//=============================================================

if(EdE==4) cout << "T = " << T*am/R << endl;
if(EdE!=4) cout << "T = " << T*am/bm/R << endl;
cout << "=======================================\n" << endl;

    //****************************START WRITING FILE WITH f function of mole fraction of density*****************
    if(p<1)
    {
        xpp_out << ";"; //Cell A1 clear
        xpu_out << ";"; //Cell A1 clear
        xpu1_out << ";"; //Cell A1 clear
        xpu2_out << ";"; //Cell A1 clear

        //xpu1_out << 0 << ";"; //Cell A2 0
        //xpu2_out << 0 << ";"; //Cell A2 0

        //Write density values on first line
        for(int i=0; i<n; i++)
        {
        xpp_out << rho_vec_out[i] << ";";//Handle P with mole fraction and density
        xpu_out << rho_vec_out[i] << ";";//Handle u with mole fraction and density
        xpu1_out << rho_vec_out[i] << ";";//Handle u with mole fraction and density
        xpu2_out << rho_vec_out[i] << ";";//Handle u with mole fraction and density
        }

        xpp_out << endl; //Jump to next line, start mole fractions P and u
        xpu_out << endl; //Jump to next line, start mole fractions P and u
        xpu1_out << endl; //Jump to next line, start mole fractions P and u
        xpu2_out << endl; //Jump to next line, start mole fractions P and u
    }

    for(i=0;i<n;i++)
    {
    xpp_out << x(0) << ";"; //Write mole fraction on first column
    xpu_out << x(0) << ";"; //Write mole fraction on first column
    xpu1_out << x(0) << ";"; //Write mole fraction on first column
    xpu2_out << x(0) << ";"; //Write mole fraction on first column DOUBT

    for(j=0; j<n; j++)
    {
        xpp_out << P_mat[i][j] << ";";//Handle P with mole fraction and density
        xpu_out << u_vec[i] << ";";//Handle u with mole fraction and density
        xpu1_out << fn1[i][j] << ";";//Handle u with mole fraction and density
        xpu2_out << fn2[i][j] << ";";//Handle u with mole fraction and density
    }

        xpp_out << bm << endl; //Jump to next line, next mole fraction
        xpu_out << bm << endl; //Jump to next line, next mole fraction
        xpu1_out << b(0) << endl; //Jump to next line, next mole fraction
        xpu2_out << b(1) << endl; //Jump to next line, next mole fraction
    }
    //****************************END WRITING FILE WITH f function of mole fraction of density*****************

k++;
p++;
    //if(x(0)==0.000001) x(0)=x(0)-0.000001;
    //} //end for x loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    T = T + step;
    g++;
} //end while T loop

    if(cp[0] != cp[1])
    {
        d2P = new double *[n];
        for(int k = 0; k <n; k++)
            d2P[k] = new double[n];

        d2u1 = new double *[n];
        for(int k = 0; k <n; k++)
            d2u1[k] = new double[n];

        d2u2 = new double *[n];
        for(int k = 0; k <n; k++)
            d2u2[k] = new double[n];

        Pmat = new double *[n];
        for(int k = 0; k <n; k++)
            Pmat[k] = new double[n];

        umat = new double *[n];
        for(int k = 0; k <n; k++)
            umat[k] = new double[n];

        u1mat = new double *[n];
        for(int k = 0; k <n; k++)
            u1mat[k] = new double[n];

        u2mat = new double *[n];
        for(int k = 0; k <n; k++)
            u2mat[k] = new double[n];

        rho_rv.resize(n);

        d2Pgen(d2P,n);
        cout << "d2P ok" << endl;
        d2u1gen(d2u1,n);
        d2u2gen(d2u2,n);
        cout << "d2u1 and d2u2 ok" << endl;
        renorm_mat_reader(Pmat,umat,n);
        cout << "Pmat umat ok" << endl;
        renorm_uu_reader(u1mat,u2mat,n);
        cout << "u1mat u2mat ok" << endl;
        rho_rv = renorm_rhovec(n);
        cout << "rho_rv ok" << endl;
        x_rv = renorm_xvec(n);
        cout << "x_rv ok" << endl;
        T = T_original;
        Renormalization=2;
        if(EdE==1) EdE = 5;
        if(EdE==3) EdE = 6;
    }
    else envelope_tracer(1e-5,env_type,n);
    break;

    //***************************************************************************
    //
    //                  CASE 3
    //                  FIND TC FAST
    //                  FIND TC FAST
    //                  FIND TC FAST
    //                  CASE 3
    //
    //****************************************************************************//



    case 3: //Find critical point using inflexion
    {
        int estim_L;
        double step_L, final_L, T_original2, Tnew_original, x_iter, T_L;
        double Tcl, Pcl, ucl, rhocl, ucline, betal, betal2, deltal;
        cout << "range L parameter? \n 1.Yes\n 2.No \n If yes, give step value after answer, and then, final L" << endl;
        cin >> estim_L;
        if(estim_L==1)
        {
        cout << "step L: ";
        cin >> step_L;
        cout << "final L: ";
        cin >> final_L;
        }
        else
        {
        final_L = L;
        step_L = L;
        }
        T_original2 = T;
        Tnew_original = final_T;

        ofstream critical_line("../critical_line.csv");
        critical_line << "Critical Line Report-----------------------------" << endl;
        critical_line << "Equation of state = " << EdE << endl;
        critical_line << "Component 1 = " << cp[0] << endl;
        critical_line << "Component 2 = " << cp[1] << endl;
        critical_line << "x1" << ";" << "Tc" << ";" << "Pc" << ";" << "rhoc" << ";" << "uc" << ";" << "beta" << ";" << "beta2" << "delta"
                      << ";" << "L" << ";" << "phi_r" << endl;

        if(cp[0]==cp[1]) x(0) = 0.999999; //If both components are equal, then consider pure component
        else
        {
            if(iter_choice==2) x(0) = 0.000;
            else
            {
            cout << "Give molar fraction for component 1: \n";
            cin >> x_iter;
            }
        }

                for(x(0)=0.00; x(0)<=1.00 ;x(0)=x(0)+0.05)
                {
                if(cp[0]==cp[1]) x(0) = 0.999999;
                cout << "BEGIN X ITERATION === x0 = " << x(0) << endl;
                if(cp[0]!=cp[1])
                {
                if(x(0)>0.000) T = Tcl - 2;
                else T = T_original2;
                }
                T_L = T;
                final_T = Tnew_original;
                if(iter_choice==1) x(0) = x_iter;

                x(1) = 1 - x(0);
                //RG parameters
                L = x(0)*pow(L_rg(0),3)+x(1)*pow(L_rg(1),3);
                L = pow(L,(1./3.));
                phi_r = x(0)*phi_rg(0) + x(1)*phi_rg(1);

        while(L<=final_L)
        {
        T = T_L;
        flag = 0;

    {
    ofstream Envelope("../Planilhas de análise/env.csv");
    ofstream envelope_exponent("../Planilhas de análise/env_exponent.csv");
    ofstream envelope_exponent2("../Planilhas de análise/env_exponent2.csv");
    ofstream crit_isot("../Planilhas de análise/crit_isotherm.csv");
    Envelope << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << ";" << "delta_P" << ";" << "delta_u" << endl;
    envelope_exponent << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << endl;
    envelope_exponent2 << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << endl;
    double Tcritical, Pcritical, rhocritical1, rhocritical2, rhocritical, ucritical;
    int counter, o;
    int exp_ok = 0;
    double Tnew;
    double b_crit, b_crit2, d_crit;
    double drho_new, drho_old;
    counter = 0;

    Tnew = T;

while(T<=Tnew)
{
    cout << "BEGIN==============================================================\n";
    //Recalculating Temperature-dependent parameters

    //DIMENSIONLESS!!!************************************************************
    //T = T_original + g*step_original;
    //final_T = final_T_original;
    //step = step_original;

    cout << "T = " << T << endl;
    if(EdE==1)
    {
        if(Tc_virtual[0] != 0)
        {
        Tr = T*Tc_virtual.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega_virtual, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, bCPA, EdE);
        }

        else
        {
        Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
        }
    }

    else
    {
    Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    if(EdE==4)
    {
        eps = 8.75;
        lambda = 1.65;
        sigma = 4.86e-10;
        L = 0.8e-4; //omega
        pi = 3.14159265359796;
        NA = 6.023e23;
        cnst = pi/6;
        l12 = 0;

        a = a_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0, T);
        b = b_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE, T);

        sig1 = pow((b(0)/cnst),(1.0/3.0));
        sig2 = pow((b(1)/cnst),(1.0/3.0));
        sig12 = 0.5*(sig1+sig2)*(1.0-l12);
        bm = (pi/6)*pow(sig12,3.0);
        am = pow((a(0)*a(1)),0.5)*(1.0-k12);

        zeta_squared = pow(lambda,2)*pow(sigma,2)/5;

        rho_max = 6/(pi*pow(sigma,3))/NA;

        //DIMENSIONLESS
        rho_max = rho_max*bm;

        phi_r = zeta_squared*10;
    }

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    w = 0;

    //DIMENSIONLESS!!!************************************************************
    T = T*bm*R/am;
    final_T = final_T*bm*R/am;
    step = step*bm*R/am;


    //====================================================================
    if(EdE==4)
    {
        T = T/bm;
        final_T = final_T/bm;
        step = step/bm;
    }

    cout << "Before renormalization =======================================\n";

if(r_type==1)
{
    //Calcular vetor de f em f0 com um cálculo
    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;


    for(k=0; k<n; k++)
    {
    rho_vec[k] = double(k)/n/bm;
    rho_vector(k) = double(k)/n/bm;
    rho_vec[0] = 1e-6;
    rho_vector(0) = 1e-6;

    //NON BONDED FRACTION DIMENSIONAL**************************
    T = T/bm/R*am;

    if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_vector(k), deltaV, X, 0, a, &Q_func, BETCR, E_auto,
                                          beta_auto, Dij);
    //cout << "X = \n" << X << endl;
    T = T*bm*R/am;
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    rho_vec[k] = rho_vec[k]*bm;
    rho_vector(k) = rho_vector(k)*bm;
    rho_vec[0] = rho_vec[0]*bm;
    rho_vector(0) = rho_vector(0)*bm;

    fv(k) = helmholtz_repulsive(EdE, R, T, rho_vector(k), am, bm, X, x, sigma, eps, kB);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << " " << x(1) << " " << am << " " << bm << endl;

    //DIMENSIONLESS!!!************************************************************
    //if(EdE != 4) fv(k) = fv(k) + 0.5*am*rho_vector(k)*rho_vector(k);
    if(EdE != 4) fv(k) = fv(k) + 0.5*rho_vector(k)*rho_vector(k);
    //DIMENSIONLESS!!!************************************************************

    f_originalv(k) = fv(k);

    //DIMENSIONLESS!!!************************************************************
    //if(EdE != 4) f_originalv(k) = fv(k) - 0.5*am*rho_vector(k)*rho_vector(k);
    if(EdE != 4) f_originalv(k) = fv(k) - 0.5*rho_vector(k)*rho_vector(k);
    //DIMENSIONLESS!!!************************************************************


    f0_vec[k] = f_originalv(k);
    f_before(k) = fv(k);
    rho = rho + rho_max/n;
    }

    //Iteração principal - o n
    for(i=1; i<9; i++) //i começando em 1
    {
        //Calcular K
        Kn = kB*T/((pow(2,3*i))*pow(L,3));

        //DIMENSIONLESS!!!************************************************************
        Kn = T/pow(2,3*i)/(pow(L,3)/bm*6.023e23);

        //DIMENSIONLESS FOR MSA!!!!!
        if(EdE==4)
        {
        T = T/R*am;

        Kn = R*T/((pow(2,3*i))*L);
        Kn = Kn/am*bm;

        T = T*R/am;
        }

        //Preparar vetores f_l e f_s
        rho = 1e-6;

        //DIMENSIONLESS!!!************************************************************
        rho = 1e-6*bm;

        if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

        for(w=0; w<n; w++)
        {
        flv(w) = helmholtz_recursion_long(EdE, fv(w), rho_vector(w), am, bm, R, T);
        fsv(w) = helmholtz_recursion_short(EdE, fv(w), rho_vector(w), am, bm, i, L, phi_r, sr_type, R, T);
        rho = rho + rho_max/n;
        }

        //cout << "flv before 201 / 201-100 = " << flv(201) << " / " << flv(201-100) << endl;
        //cout << "fsv before 201 / 201-100 = " << fsv(201) << " / " << fsv(201-100) << endl;

            //Iteração 2 - calcular os valores para f no i atual
            rho = 1e-6;

            //DIMENSIONLESS!!!************************************************************
            rho = 1e-6*bm;

            w = 0;
            width = rho_max/n;

            for(w=0; w<n; w++)
            {
                if(EdE != 4) delta_fv(w) = df_calculation(w,n,Iteration,width,Kn,rho_vec,flv,fsv);
                if(EdE == 4)
                {
                    if(w<n/2) delta_fv(w) = df_calculation(w,n,Iteration,width,Kn,rho_vec,flv,fsv);
                }
            }

        //Calcular o novo vetor de f, ajustando com o vetor de delta_f
        fv.array() = fv.array() + delta_fv.array();
        cout << "i = " << i << " / f 201 = " << fv(201) << " / delta 201 --> " << delta_fv(201) << endl;
        //cin >> stop;
    }

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    for(w=0; w<n; w++)
    {
    //DIMENSIONLESS!!!************************************************************
    //if(EdE != 4) fv(w) = fv(w) - 0.5*am*rho_vector(w)*rho_vector(w);
    if(EdE != 4) fv(w) = fv(w) - 0.5*rho_vector(w)*rho_vector(w);
    //DIMENSIONLESS!!!************************************************************

    f_vec[w] = fv(w);

    //DIMENSIONLESS!!!************************************************************
    if(EdE != 4)
    {
    f_vec[w] = fv(w)*am/bm/bm;
    f0_vec[w] = f0_vec[w]*am/bm/bm;
    rho_vec[w] = rho_vec[w]/bm;
    }

    if(EdE == 4)
    {
    f_vec[w] = fv(w)*am/bm;
    f0_vec[w] = f0_vec[w]*am/bm;
    rho_vec[w] = rho_vec[w]/bm;
    }

    //cout << "rho = " << rho_vector(w) << "  //  f = " << fv(w) << endl;
    Not_splined << std::fixed << std::setprecision(15) << rho_vector(w) << ";" << fv(w) << ";" << f_originalv(w) << ";"
                << T << endl;
    rho = rho + rho_max/n;
    }
}
    //====================================================================

rho = 1e-6;
w=0;


for(w=0; w<n; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/n/bm;
    rho_vec_out[0] = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);


for(i=0; i<n; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b1;
    a = (P_vec[n/2+5]-P_vec[n/2-5])/(rho_vec_out[n/2+5]-rho_vec_out[n/2-5]);
    b1 = P_vec[n/2+5]-a*rho_vec_out[n/2+5];
    P_vec[n/2-4] = a*rho_vec_out[n/2-4] + b1;
    P_vec[n/2-3] = a*rho_vec_out[n/2-3] + b1;
    P_vec[n/2-2] = a*rho_vec_out[n/2-2] + b1;
    P_vec[n/2-1] = a*rho_vec_out[n/2-1] + b1;
    P_vec[n/2] = a*rho_vec_out[n/2] + b1;
    P_vec[n/2+1] = a*rho_vec_out[n/2+1] + b1;
    P_vec[n/2+2] = a*rho_vec_out[n/2+2] + b1;
    P_vec[n/2+3] = a*rho_vec_out[n/2+3] + b1;
    P_vec[n/2+4] = a*rho_vec_out[n/2+4] + b1;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<n; i++)
{

    if(EdE!=4)
    {
    Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i]<< ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/bm/R << endl;
    }

    if(EdE==4)
    {
    Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i]<< ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/R << endl;
    }

}

//=============================================================


for(i=0;i<n;i++)
    {
        P_vec_e(i) = bm*bm*P_vec[i]/am;
        u_vec_e(i) = bm*u_vec[i]/am;
        rho_vec_out_e(i) = bm*rho_vec_out[i];
        f_vec_e(i) = -P_vec_e(i) + rho_vec_out_e(i)*u_vec_e(i);

        V_vec[i] = 1/rho_vec_out[i];
        A_vec[i] = f_vec_out[i]/rho_vec_out[i];

        P_env[i] = P_vec_e(i);
        rho_env[i] = rho_vec_out_e(i);
        f_env[i] = f_vec_e(i);
        u_env[i] = u_vec_e(i);
    }

        std::vector<double> dta(6), out(2);
        T = T/bm/R*am;
        if(EdE==4) T = T*bm;
        cout << "T = " << T << endl;
        cout << "=======================================\n" << endl;

        k++;
        g++;

        std::vector<double> d1P(n);
        d1P = fin_diff_1_vec(rho_vec_out, P_vec);
        int inflex;
        for(w=1; w<n-10; w++)
        {
            if(d1P[w] < 0)
            {
                inflex = 1; //Isotherm has inflexion point
                w = n;
            }

            else inflex = 0; //Isotherm does not have inflexion point
        }

        if(flag==6)
        {
            dta = dens_newt(rho_vec_out,f_vec_out,P_vec,u_vec,1e-5,n);
            cout << "\n" << dta[0] << " / " << dta[1] << " / " << dta[2] << " / " << dta[3] << " / " << dta[4] << endl;

            Tcritical = T;
            Pcritical = dta[2];
            rhocritical1 = dta[0];
            rhocritical2 = dta[1];
            rhocritical = (rhocritical1+rhocritical2)/2;
            ucritical = dta[4];
            flag = 7;
        }

        if(inflex==1 && flag<5)
        {
            switch(flag)
            {
                case 0: T = T + 1;
                break;

                case 1: T = T + 0.1;
                break;

                case 2: T = T + 0.01;
                break;

                case 3: T = T + 0.001;
            }
            Tnew = T;
        }

        if(inflex==0 && flag<5)
        {
            out = T_tracer_inflexion(EdE, T, flag);
            T = out[0];
            flag = out[1];
        }

        if(flag==9)
        {
            dta = dens_newt(rho_vec_out,f_vec_out,P_vec,u_vec,1e-5,n);
            cout << "\n" << dta[0] << " / " << dta[1] << " / " << dta[2] << " / " << dta[3] << " / " << dta[4] << endl;
            if(dta[0] > 0) envelope_exponent2 << T << ";" << dta[0] << ";" << dta[1] << ";" << dta[2] << ";" << dta[4] << endl;

            if(T>=0.998*Tcritical)
            {
                cout << "will calculate beta2!" << endl;
                exp_ok = 2;
                T = Tnew*10;
            }

            cout << "\nflag 9 | " << " T = " << T << " / Tc = " << Tcritical << " \ " << exponents << "\ " << exp_ok << endl;

            if(exponents==1 && exp_ok==2)
            {
            //double beta_crit = critical_exponents(EdE);
            beta_exponent2(&b_crit2, Tcritical, rhocritical);
            cout << "beta2 out of critical = " << b_crit2 << endl;
            T = Tnew+1;
            }


            T = T + 0.001*Tcritical;
            cout << "fnew T to calculate beta 2 = " << T << endl;
            o++;
        }

        if(flag==8)
        {

            dta = dens_newt(rho_vec_out,f_vec_out,P_vec,u_vec,1e-5,n);
            cout << "\n" << dta[0] << " / " << dta[1] << " / " << dta[2] << " / " << dta[3] << " / " << dta[4] << endl;
            if(dta[0] > 0) envelope_exponent << T << ";" << dta[0] << ";" << dta[1] << ";" << dta[2] << ";" << dta[4] << endl;


            if(T>=Tcritical)
            {
                cout << "will calculate beta!" << endl;
                exp_ok = 1;
                T = Tnew*10;
            }

            cout << "\nflag 8 | " << " T = " << T << " / Tc = " << Tcritical << " \ " << exponents << "\ " << exp_ok << endl;

            if(exponents==1 && exp_ok==1)
            {
            //double beta_crit = critical_exponents(EdE);
            beta_exponent(&b_crit,Tcritical,rhocritical);
            cout << "beta out of critical = " << b_crit << endl;
            T = Tnew+1;
            flag = 9;
            Tnew = 1e5;
            }

            T = T + 0.1;

            if(flag == 9)
            {
                T = T - 0.1;
                T = 0.985*Tcritical;
                cout << "first T to calculate beta 2 = " << T << endl;
            }
        }

        if(flag==7)
        {

        for(int l=0; l<n; l++)
        {
        crit_isot << std::fixed << std::setprecision(15) << rho_vec_out[l] << ";" << f_vec_out[l] << ";"
                  << f0_vec_out[l] << ";" << u_vec[l] << ";" << P_vec[l]<< ";" << u_vec_0[l] << ";" << P_vec_0[l] << ";" << T << endl;
        }

        cout << "\nflag 7" << endl;
        flag = 8;

            Tcritical = T;
            Pcritical = dta[2];
            rhocritical1 = dta[0];
            rhocritical2 = dta[1];
            rhocritical = (rhocritical1+rhocritical2)/2;
            ucritical = dta[4];

        delta_exponent(&d_crit, rhocritical, ucritical, n);
        cout << "delta out of critical = " << d_crit << endl;
        cin >> g;

        T = T - 1;
        Tnew = 2*T;
        }


}
         ofstream out_simulation("record.csv", fstream::app);
         out_simulation << EdE << ";" << cp[0] << ";" << Tcritical << ";" << Pcritical << ";" << rhocritical << ";"
                        << ucritical << ";" << b_crit << ";" << b_crit2 << ";" << d_crit << ";" << L << ";" << phi_r << endl;

        Tcl = Tcritical;
        Pcl = Pcritical;
        rhocl = rhocritical;
        ucline = ucritical;
        betal = b_crit;
        betal2 = b_crit2;
        deltal = d_crit;

        Envelope.close();
        envelope_exponent.close();
        envelope_exponent2.close();
        crit_isot.close();

        if(L + step_L <= final_L)
        {
        //Clean Data
        Envelope.open("../Planilhas de análise/env.csv", std::fstream::out | std::fstream::trunc);
        envelope_exponent.open("../Planilhas de análise/env_exponent.csv", std::fstream::out | std::fstream::trunc);
        envelope_exponent2.open("../Planilhas de análise/env_exponent2.csv", std::fstream::out | std::fstream::trunc);
        crit_isot.open("../Planilhas de análise/crit_isotherm.csv", std::fstream::out | std::fstream::trunc);

        //Close again
        Envelope.close();
        envelope_exponent.close();
        envelope_exponent2.close();
        crit_isot.close();
        }

    }

        L = L + step_L;
        }

                critical_line << Tcl << ";" << Pcl << ";" << rhocl << ";"
                              << ucline << ";" << betal << ";" << betal2 << ";" << deltal
                              << ";" << L-step_L << ";" << phi_r << endl;

                if(iter_choice==1) x(0) = 1; //Stop iteration, user chose to calculate single point
                }

    }
    break;

    case 4: //Binary mixture
{
    double **bfnl;
    bfnl = new double *[n];
    for(int k = 0; k <n; k++)
        bfnl[k] = new double[n];

    double **bfns;
    bfns = new double *[n];
    for(int k = 0; k <n; k++)
        bfns[k] = new double[n];

    double **rho_mat;
    rho_mat = new double *[n];
    for(int k = 0; k <n; k++)
        rho_mat[k] = new double[n];

    double **fv_mat;
    fv_mat = new double *[n];
    for(int k = 0; k <n; k++)
        fv_mat[k] = new double[n];

    double **f0_mat;
    f0_mat = new double *[n];
    for(int k = 0; k <n; k++)
        f0_mat[k] = new double[n];

    double **f_originalm;
    f_originalm = new double *[n];
    for(int k = 0; k <n; k++)
        f_originalm[k] = new double[n];

    double **f_mat;
    f_mat = new double *[n];
    for(int k = 0; k <n; k++)
        f_mat[k] = new double[n];

    double **flm;
    flm = new double *[n];
    for(int k = 0; k <n; k++)
        flm[k] = new double[n];

    double **fsm;
    fsm = new double *[n];
    for(int k = 0; k <n; k++)
        fsm[k] = new double[n];

    double di, dfm, dy;
    double rho1a[n], rho2a[n];

    dy = double(1/n);
    cout << "dy inicial = " << dy << " " << 1/n << " " << double(1/n) << endl;
    dy = pow(n,-1);
    cout << "dy inicial = " << dy << " " << 1/n << " " << double(1/n) << endl;

    while(T<=final_T)
    {
    int p=0;

    for(x(0)=0.000; x(0)<=1.001; x(0) = x(0) + 0.005)
    {
    cout << "x0 = " << x(0) << endl;

    if(cp[0]==cp[1]) x(0) = 1; //If both components are equal, then consider pure component

    x(1) = 1 - x(0);
    //RG parameters
    L = x(0)*pow(L_rg(0),3)+x(1)*pow(L_rg(1),3);
    L = pow(L,(1./3.));
    phi_r = x(0)*phi_rg(0) + x(1)*phi_rg(1);
    cout << "BEGIN==============================================================\n";
    //Recalculating Temperature-dependent parameters

    //DIMENSIONLESS!!!************************************************************
    T = T_original + g*step_original;
    final_T = final_T_original;
    step = step_original;

    cout << "T = " << T << endl;
    if(EdE==1)
    {
        if(Tc_virtual[0] != 0)
        {
        Tr = T*Tc_virtual.asDiagonal().inverse().diagonal(); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega_virtual, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega_virtual, Tc_virtual, Pc_virtual, alfa, bCPA, EdE);
        }

        else
        {
        Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
        alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
        a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
        b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
        }
    }

    else
    {
    Tr = T*(Tc.asDiagonal().inverse().diagonal()); //Vetor com temperaturas reduzidas virtuais
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    cout << "n = " << n << endl;

    if(EdE==4)
    {
        eps = 8.75;
        lambda = 1.65;
        sigma = 4.86e-10;
        L = 0.8e-4; //omega
        pi = 3.14159265359796;
        NA = 6.023e23;
        cnst = pi/6;
        l12 = 0;

        a = a_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0, T);
        b = b_function_msa(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE, T);

        sig1 = pow((b(0)/cnst),(1.0/3.0));
        sig2 = pow((b(1)/cnst),(1.0/3.0));
        sig12 = 0.5*(sig1+sig2)*(1.0-l12);
        bm = (pi/6)*pow(sig12,3.0);
        am = pow((a(0)*a(1)),0.5)*(1.0-k12);
        /*
        cout << "sig1 = " << sig1 << endl;
        cout << "sig2 = " << sig2 << endl;
        cout << "sig12 = " << sig12 << endl;
        cout << "b(0) = " << b(0) << endl;
        cout << "b(1) = " << b(1) << endl;
        cout << "cnst = " << cnst << endl;
        cout << "bm = " << bm << endl;
        cout << "am = " << am << endl;
        */

        zeta_squared = pow(lambda,2)*pow(sigma,2)/5;

        rho_max = 6/(pi*pow(sigma,3))/NA;

        //DIMENSIONLESS
        rho_max = rho_max*bm;

        phi_r = zeta_squared*10;
    }

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    w = 0;

    //DIMENSIONLESS!!!************************************************************
    T = T*bm*R/am;
    final_T = final_T_original*bm*R/am;
    step = step_original*bm*R/am;

    if(EdE==4)
    {
        T = T/bm;
        final_T = final_T/bm;
        step = step/bm;
    }


    //====================================================================
    cout << "Before renormalization =======================================\n";

if(r_type==1)
{
    //Calcular vetor de f em f0 com um cálculo
    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    for(int k1=0; k1<n; k1++)
    {
    rho1_vec[k1] = double(k1)/n/b(0);
    rho2_vec[k1] = double(k1)/n/b(0);
    rho1_vec[0] = 1e-6;
    rho2_vec[0] = 1e-6;

    rho1_vec[k1] = rho1_vec[k1]*b(0);
    rho2_vec[k1] = rho2_vec[k1]*b(1);

    cout << "rho 1 / 2 = " << rho1_vec[k1] << " " << rho2_vec[k1] << " " << k1 << endl;
    }

    for(int k1=0; k1<n; k1++)
    {
        for(int k2=0; k2<n; k2++)
        {
        x(0) = double(k1)/n;
        x(1) = double(k2)/n;


    //NON BONDED FRACTION DIMENSIONAL**************************
    T = T/bm/R*am;

    if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_vector(k), deltaV, X, 0, a, &Q_func, BETCR, E_auto,
                                          beta_auto, Dij);

    T = T*bm*R/am;
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    rho_mat[k1][k2] = x(0)*rho_vec[k1] + x(1)*rho2_vec[k2];

    fv_mat[k1][k2] = helmholtz_repulsive(5, R, T, rho_mat[k1][k2], am, bm, X, x, sigma, eps, kB);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << endl;

    f_originalm[k1][k2] = fv_mat[k1][k2];

    if(EdE != 4) fv_mat[k1][k2] = fv_mat[k1][k2] + 0.5*rho_mat[k1][k2]*rho_mat[k1][k2];

    f0_mat[k1][k2] = f_originalm[k1][k2];
    rho = rho + rho_max/n;
        }
    }
    cout << "fv_mat ok" << endl;

    //Iteração principal - o n
    for(i=1; i<9; i++) //i começando em 1
    {
        cout << "renormalizing: i = " << i << endl;
        //Calcular K
        Kn = kB*T/((pow(2,3*i))*pow(L,3));

        //DIMENSIONLESS!!!************************************************************
        Kn = T/pow(2,3*i)/(pow(L,3)/bm*6.023e23);


        //DIMENSIONLESS FOR MSA!!!!!
        if(EdE==4)
        {
        T = T/R*am;

        Kn = R*T/((pow(2,3*i))*L);
        Kn = Kn/am*bm;

        T = T*R/am;
        }

        //Kn = kB*T*NA/((pow(2,3*i))*L*b(0));

        //Preparar vetores f_l e f_s
        rho = 1e-6;

        //DIMENSIONLESS!!!************************************************************
        rho = 1e-6*bm;

        if(EdE==3 || EdE==6) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto, Dij);

        for(int k1=0; k1<n; k1++)
        {
            for(int k2=0; k2<n; k2++)
            {
            flm[k1][k2] = helmholtz_recursion_long(EdE, fv_mat[k1][k2], rho_mat[k1][k2], am, bm, R, T);
            fsm[k1][k2] = helmholtz_recursion_short(EdE, fv_mat[k1][k2], rho_mat[k1][k2], am, bm, i, L, phi_r, sr_type, R, T);
            rho = rho + rho_max/n;
            bfnl[k1][k2] = flm[k1][k2];
            bfns[k1][k2] = fsm[k1][k2];
            }
        }

        cout << "calculated bfnl and bfns" << endl;

        //cout << "flv before 201 / 201-100 = " << flv(201) << " / " << flv(201-100) << endl;
        //cout << "fsv before 201 / 201-100 = " << fsv(201) << " / " << fsv(201-100) << endl;

            //Iteração 2 - calcular os valores para f no i atual
            rho = 1e-6;

            //DIMENSIONLESS!!!************************************************************
            rho = 1e-6*bm;

            w = 0;
            width = rho_max/n;

            for(int k1=0; k1<n; k1++)
            {
                for(int k2=0; k2<n; k2++)
                {
                if(EdE != 4) di = di_calculation(bfnl,bfns,1/Kn,dy,k1,k2,n);
                //cout << "di = " << di << endl;
                    if(EdE == 4)
                    {
                    if(rho_mat[k1][k2]<n/2) di = di_calculation(bfnl,bfns,1/Kn,dy,k1,k2,n);
                    else di = 0;
                    }
                dfm = -di*Kn;
                cout << "k1,k2: " << k1 << " " << k2 << " " << di << " " << dfm << " " << fv_mat[k1][k2] << " " << Kn << endl;
                fv_mat[k1][k2] = fv_mat[k1][k2] + dfm;
                }
            }

            cout << "corrected fv_mat" << endl;


        //Calcular o novo vetor de f, ajustando com o vetor de delta_f
        //fv.array() = fv.array() + delta_fv.array();
        cout << "i = " << i << " / Kn = " << Kn*am/bm << " / f 200/200 = " << fv_mat[200][200]*am/bm << " / delta 60 --> "
             << delta_fv(60)*am/bm << endl;
        //cin >> stop;
        /*
        for(w=0; w<n; w++)
        {
            dfnout << i << ";" << rho_vector(w)/bm << ";" << fv(w)*am/bm << ";" << delta_fv(w)*am/bm << endl;
        }
        */
    }

    rho = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    rho = 1e-6*bm;

    for(int k1=0; k1<n; k1++)
    {
        for(int k2=0; k2<n; k2++)
        {
    if(EdE != 4) fv_mat[k1][k2] = fv_mat[k1][k2] - 0.5*rho_mat[k1][k2]*rho_mat[k1][k2];

    f_mat[k1][k2] = fv_mat[k1][k2];

    //DIMENSIONLESS!!!************************************************************
    if(EdE != 4)
    {
    f_vec[w] = fv(w)*am/bm/bm;
    f0_vec[w] = f0_vec[w]*am/bm/bm;
    rho_vec[w] = rho_vec[w]/bm;
    }

    if(EdE == 4)
    {
    f_vec[w] = fv(w)*am/bm;
    f0_vec[w] = f0_vec[w]*am/bm;
    rho_vec[w] = rho_vec[w]/bm;
    }


    //cout << "rho = " << rho_vector(w) << "  //  f = " << fv(w) << endl;
    Not_splined << std::fixed << std::setprecision(15) << rho_vector(w) << ";" << fv(w) << ";" << f_originalv(w) << ";"
                << T << endl;
    rho = rho + rho_max/n;
        }
    }
}
    //====================================================================

rho = 1e-6;
w=0;

for(w=0; w<n; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/n/bm;
    rho_vec_out[0] = 1e-6;

    rho1_vec[w] = double(w)/n/b(0);
    rho2_vec[w] = double(w)/n/b(1);
    rho1_vec_out[w] = double(w)/n/b(0);
    rho2_vec_out[w] = double(w)/n/b(1);

    rho1a[w] = rho1_vec[w];
    rho2a[w] = rho2_vec[w];

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

    if(EdE!=4) T = T*am/bm/R;
    if(EdE==4) T = T*am/R;

//Subtract ideal gas contribution before cubic spline
for(w=0; w<n; w++)
{
    f_vec[w] = f_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f0_vec[w] = f0_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
}

f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);

u_vec_res1 = cspline_deriv1_vec(rho1_vec_out, f_vec, rho1_vec_out);
u_vec_res2 = cspline_deriv1_vec(rho2_vec_out, f_vec, rho2_vec_out);

//Add ideal gas contribution after cubic spline
for(w=0; w<n; w++)
{
    f_vec_out[w] = f_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    f0_vec_out[w] = f0_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    u_vec[w] = u_vec[w] + R*T*(log(rho_vec_out[w]));
    u_vec_0[w] = u_vec_0[w] + R*T*(log(rho_vec_out[w]));
}
    if(EdE!=4) T = T/am*bm*R;
    if(EdE==4) T = T/am*R;

for(i=0; i<n; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b1;
    a = (P_vec[n/2+5]-P_vec[n/2-5])/(rho_vec_out[n/2+5]-rho_vec_out[n/2-5]);
    b1 = P_vec[n/2+5]-a*rho_vec_out[n/2+5];
    P_vec[n/2-4] = a*rho_vec_out[n/2-4] + b1;
    P_vec[n/2-3] = a*rho_vec_out[n/2-3] + b1;
    P_vec[n/2-2] = a*rho_vec_out[n/2-2] + b1;
    P_vec[n/2-1] = a*rho_vec_out[n/2-1] + b1;
    P_vec[n/2] = a*rho_vec_out[n/2] + b1;
    P_vec[n/2+1] = a*rho_vec_out[n/2+1] + b1;
    P_vec[n/2+2] = a*rho_vec_out[n/2+2] + b1;
    P_vec[n/2+3] = a*rho_vec_out[n/2+3] + b1;
    P_vec[n/2+4] = a*rho_vec_out[n/2+4] + b1;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<n; i++)
{
// Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
//        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;

        //DIMENSIONLESS!!!************************************************************
 if(EdE!=4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/bm/R << endl;
 }

 if(EdE==4)
 {
 Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i]<< ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T*am/R << endl;
 }

}

//=============================================================

if(EdE==4) cout << "T = " << T*am/R << endl;
if(EdE!=4) cout << "T = " << T*am/bm/R << endl;
cout << "=======================================\n" << endl;

    //****************************START WRITING FILE WITH f function of mole fraction of density*****************
    if(p<1)
    {
        xpp_out << ";"; //Cell A1 clear
        xpu_out << ";"; //Cell A1 clear
        xpu1_out << ";"; //Cell A1 clear
        xpu2_out << ";"; //Cell A1 clear

        //Write density values on first line
        for(int i=0; i<n; i++)
        {
        xpp_out << rho_vec_out[i]*bm << ";";//Handle P with mole fraction and density
        xpu_out << rho_vec_out[i]*bm << ";";//Handle u with mole fraction and density
        xpu1_out << rho1_vec_out[i]*b(0) << ";";//Handle u with mole fraction and density
        xpu2_out << rho2_vec_out[i]*b(1) << ";";//Handle u with mole fraction and density
        }

        xpp_out << endl; //Jump to next line, start mole fractions P and u
        xpu_out << endl; //Jump to next line, start mole fractions P and u
        xpu1_out << endl; //Jump to next line, start mole fractions P and u
        xpu2_out << endl; //Jump to next line, start mole fractions P and u
    }


    xpp_out << x(0) << ";"; //Write mole fraction on first column
    xpu_out << x(0) << ";"; //Write mole fraction on first column
    xpu1_out << x(0) << ";"; //Write mole fraction on first column
    xpu2_out << x(0) << ";"; //Write mole fraction on first column DOUBT

    for(i=0; i<n; i++)
    {
        xpp_out << P_vec[i] << ";";//Handle P with mole fraction and density
        xpu_out << u_vec[i] << ";";//Handle u with mole fraction and density
        xpu1_out << u_vec_res1[i] << ";";//Handle u with mole fraction and density
        xpu2_out << u_vec_res2[i] << ";";//Handle u with mole fraction and density
    }

        xpp_out << bm << endl; //Jump to next line, next mole fraction
        xpu_out << bm << endl; //Jump to next line, next mole fraction
        xpu1_out << b(0) << endl; //Jump to next line, next mole fraction
        xpu2_out << b(1) << endl; //Jump to next line, next mole fraction
    //****************************END WRITING FILE WITH f function of mole fraction of density*****************

k++;
p++;

    } //end for x loop

    T = T + step;
    g++;
} //end while T loop

    if(cp[0] != cp[1])
    {
        d2P = new double *[200];
        for(int k = 0; k <200; k++)
            d2P[k] = new double[n];

        d2u = new double *[200];
        for(int k = 0; k <200; k++)
            d2u[k] = new double[n];

        Pmat = new double *[200];
        for(int k = 0; k <200; k++)
            Pmat[k] = new double[n];

        umat = new double *[200];
        for(int k = 0; k <200; k++)
            umat[k] = new double[n];

        u1mat = new double *[200];
        for(int k = 0; k <200; k++)
            u1mat[k] = new double[n];

        u2mat = new double *[200];
        for(int k = 0; k <200; k++)
            u2mat[k] = new double[n];

        p2 = new double *[n];
        for(int k = 0; k <n; k++)
            p2[k] = new double[n];

        d2u1 = new double *[n];
        for(int k = 0; k <n; k++)
            d2u1[k] = new double[n];

        d2u2 = new double *[n];
        for(int k = 0; k <n; k++)
            d2u2[k] = new double[n];

        rho_rv.resize(n);

        d2Pgen(d2P,n);
        d2u1gen(d2u1,n);
        d2u2gen(d2u2,n);
        renorm_mat_reader(Pmat,umat,n);
        renorm_uu_reader(u1mat,u2mat,n);
        rho_rv = renorm_rhovec(n);
        x_rv = renorm_xvec(n);

        d2_chem_p(f_mat,rho1a,rho2a,n,b(0)/n,b(1)/n,p2,d2u1,d2u1);

        T = T_original;
        Renormalization=2;
        EdE = 5;
    }
    else envelope_tracer(1e-5,env_type,n);
}
    break;


}

}

if(Renormalization==2)
{

//CÁLCULO PARA COMPONENTE PURO
//==============================================================================================================
if(mixture==1)
{

if(process==1)
{
    cout << "\nDefine initial Temperature: ";
    cin >> init_T;

    cout << "\nDefine final Temperature: ";
    cin >> final_T;

    cout << "\nDefine steps: ";
    cin >> step;

    T = init_T;
    Told = T;
}

if(process==2)
{
    cout << "\nDefine initial Pressure: ";
    cin >> init_P;

    cout << "\nDefine final Pressure: ";
    cin >> final_P;

    cout << "\nDefine steps: ";
    cin >> step;

    P = init_P;
    Pold = P;
}

if(process==1)
{
for (T = init_T ; T<=final_T ; T=T+step)
{
 x(0) = 0.999999;
 x(1) = 1-x(0);

Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);

E_row = ((E_row.array().log())*Told/T).exp();
E_col = ((E_col.array().log())*Told/T).exp();
E_auto = ((E_auto.array().log())*Told/T).exp();

Told = T;

 if(iter_choice==1)
 {
 cout << "Define Temperature: ";
 cin >> T;

 counter = 0;

    if(counter==0)
    {
        switch(process)
        {
            case 1: //Isothermic
            P = Psat.transpose()*x;
            Pinit = P;
            break;

            case 2: //Isobaric
            T = Tsat.transpose()*x;
            Tinit = T;
            CT = C.array()+T;
            logPsat = A - (CT.asDiagonal().inverse()*B);
            ln10.fill(log(10));
            lnPsat = (logPsat*ln10.transpose()).diagonal();
            Psat = lnPsat.array().exp();
            cout << "Tinit" << Tinit << endl;
            break;
            //cin.get();
        }

    y = ((Psat*x.transpose()).diagonal()).array()/P;
    yinit = y;
    }

 }

 if(nc>2)
 {
     iter_choice=1;

     int q;
     q = 0;
     while(q<nc)
     {
     cout << "Define x for component " << q+1 << " : ";
     cin >> x(q);
     q++;
     }
 }

switch(process)
{
case 1: //Isothermic
P = Pinit;

CT = C.array()+T;
//logPsat = A - (B*(CT.asDiagonal().inverse().diagonal().transpose())).diagonal();
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();
P = Psat.transpose()*x;
break;

case 2: //Isobaric
T = Tinit;
CT = C.array()+T;
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();

Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
break;
}

if(counter == max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    P = Psat.transpose()*x;
    break;

    case 2: //Isobaric
    T = Tsat.transpose()*x;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    break;
    }

y = ((Psat*x.transpose()).diagonal()).array()/P;
}

if(counter != max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    //P = Pinit;

    CT = C.array()+T;
//logPsat = A - (B*(CT.asDiagonal().inverse().diagonal().transpose())).diagonal();
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();
P = Psat.transpose()*x;

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }
    break;

    case 2: //Isobaric
    T = Tinit;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }

    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    break;
    }

Vlinit = Vl;
Vvinit = Vv;
}

counter = 0;
int k;
k=1;
errorKx = tolKx + 1;
tol_u = 0.00001;

while(errorKx>tolKx)
{
Tr = T*(Tc.asDiagonal().inverse().diagonal());
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);

double Bcpa;
MatrixXd pre_F(nc,4);
    VectorXd one_4(4), one_4nc(4*nc);
    Bcpa = x.transpose()*b;
one_4 <<    1,
            1,
            1,
            1;

for(i=0;i<(4*nc);i++)
{
    one_4nc(i) = 1;
}

    //Liquid phase fugacity calculation
    //am and bm calculation
    phase = 1; //1 for liquid, 2 for vapor
    bl = b_mixing_rules_function(nc, b, x, MR);
    al = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    phase = 2;
    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    if(counter == 0 || counter == max_num_iter)
    {
        Vlinit = bl/0.99; //iota == 0.99
        Vvinit = bv+(R*T/P); //iota = bv/(bv+(R*T/P), Vvinit = bv/iota

        if(iter_choice==1)
        {
        Vl = Vlinit;
        Vv = Vvinit;
        }
    }


    deltaV = 0;
    phase = 1;
    Xl = volume_function(nc, EdE, phase, x, Xl, EdE_parameters, bl, al, R, T, P, tolV, tolZl, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto, beta_auto, Dij);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto, Dij);

    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);

    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vl1, Vl2);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
Kx = (x.asDiagonal())*K;

    for(i=0; i<nc; i++)
    {
         one(i) = 1;
    }

sumKx = one.transpose()*Kx;
double sumKxold;
sumKxold = sumKx;
errorSUMKx = tolSUMKx + 1;


double counter2;
counter2 = 0;

while(errorSUMKx>tolSUMKx || counter2<=1)
    {
    y = Kx.array()/sumKx;

    initialSUMKx = sumKx;

    //Vapor phase fugacity calculation
    //am and bm calculation
    if(process==2)
    {
    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    phase = 2; //1 for liquid, 2 for vapor


    phase = 2;
    deltaV = 0;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto,
                        beta_auto, Dij);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);


    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


    K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
    Kx = (x.asDiagonal())*K;
    sumKxnew = one.transpose()*Kx;

    finalSUMKx = sumKxnew;
    errorSUMKx = fabs(finalSUMKx - initialSUMKx);

    errorSUMKx = errorSUMKx/finalSUMKx;

    sumKx = sumKxnew;

 if(counter2==200)
 {
   sumKx = sumKxold;
   errorSUMKx = 0.00000000000001;
 }
 counter2++;

    }

double errorKxnew;
Ey = sumKx-1;
errorKx = fabs(Ey);

errorKx = errorKx/sumKx;

y = Kx.array()/sumKx;


switch(process)
{
case 1: //Isothermic
P = P*sumKx;
break;

case 2: //Isobaric

    T = 0.1*T/sumKx+0.9*T; //AQUI DIVIDE

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

break;
}

if(isnan(errorKx)==1 && process==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
}

counter++;

double trivial, V_check;

switch(process)
{
case 1: //Isothermic
    if(isnan(P)==1 || isinf(P)==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "P = NAN OR INF" << endl;
    cin.get();
}
break;

case 2: //Isobaric
    if(isnan(T)==1 || isinf(T)==1)
{
    T = (1+0.1*k)*T;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "T = NAN OR INF" << endl;

    counter = 500;
    //cin.get();
}
break;
}

if(counter==max_num_iter)
    {
    errorKx=0.00000000000001;
    }
}

if(counter!=max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    Pinit = P;
    break;

    case 2: //Isobaric
    Tinit = T;
    break;
    }
    Vlinit = Vl;
    Vvinit = Vv;
    yinit = y;
}

if(counter==max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
        Pinit = Psat.transpose()*x;
        break;

    case 2: //Isobaric
        Tinit = Tsat.transpose()*x;
        break;
    }

    yinit = ((Psat*x.transpose()).diagonal()).array()/Pinit;
}

//------------------------------------------
//P = P*100; //Converting from bar para kPa
//Converting directly on output
y = Kx;
gama = Psat.asDiagonal().inverse()*(x.asDiagonal().inverse()*y);
gama = P*gama.array();
ln_gama = gama.array().log();
G_ex = x.transpose()*ln_gama;
G_ex = G_ex*R*T;
cout << "--------------------------------" << endl;
cout << "Zl = " << Zl << endl;
cout << "Zv = " << Zv << endl;
cout << "Vl = " << Vl << endl;
cout << "Vv = " << Vv << endl;
cout << "x1 = " << x(0) << endl;
cout << "y1 = " << y(0) << endl;
cout << "P(bar) = " << P << endl;
cout << "T(K) = " << T << endl;
cout << "errorKx = " << errorKx << endl;
cout <<"--------- counter = " << counter << " ---------" << endl;
    if(process==1)
    {
    output << x(0) << ";" << y(0) << ";" << T << ";" << P*100 << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex
           << ";" << P << ";" << 1/Vl << ";" << 1/Vv << endl;
    }

    if(process==2)
    {
    output << x(0) << ";" << y(0) << ";" << P*100 << ";" << T << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex << endl;
    }


if(iter_choice==1)
 {
 cout << "End of calculation \n \n";
 counter = 0;
 }

}

}


if(process==2)
{
for (P = init_P ; P<=final_P ; P=P+step)
{
 x(0) = 0.999999;
 x(1) = 1-x(0);

 if(iter_choice==1)
 {
 cout << "Define Pressure: ";
 cin >> P;

 counter = 0;

    if(counter==0)
    {
        switch(process)
        {
            case 1: //Isothermic
            P = Psat.transpose()*x;
            Pinit = P;
            break;

            case 2: //Isobaric
            T = Tsat.transpose()*x;
            Tinit = T;
            CT = C.array()+T;
            logPsat = A - (CT.asDiagonal().inverse()*B);
            ln10.fill(log(10));
            lnPsat = (logPsat*ln10.transpose()).diagonal();
            Psat = lnPsat.array().exp();
            cout << "Tinit" << Tinit << endl;
            break;
            //cin.get();
        }

    y = ((Psat*x.transpose()).diagonal()).array()/P;
    yinit = y;
    }

 }

 if(nc>2)
 {
     iter_choice=1;

     int q;
     q = 0;
     while(q<nc)
     {
     cout << "Define x for component " << q+1 << " : ";
     cin >> x(q);
     q++;
     }
 }

/*
switch(process)
{
case 1: //Isothermic
P = Pinit;
break;

case 2: //Isobaric
T = Tinit;
CT = C.array()+T;
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();

Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
break;
}
*/

if(counter == max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    P = Psat.transpose()*x;
    break;

    case 2: //Isobaric
    T = Tsat.transpose()*x;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    break;
    }

y = ((Psat*x.transpose()).diagonal()).array()/P;
}

if(counter != max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    P = Pinit;

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }
    break;

    case 2: //Isobaric
    T = Tinit;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }

    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    break;
    }

Vlinit = Vl;
Vvinit = Vv;
}

counter = 0;
int k;
k=1;
errorKx = tolKx + 1;
tol_u = 0.00001;

while(errorKx>tolKx)
{
Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);

double Bcpa;
MatrixXd pre_F(nc,4);
    VectorXd one_4(4), one_4nc(4*nc);
    Bcpa = x.transpose()*b;
one_4 <<    1,
            1,
            1,
            1;

for(i=0;i<(4*nc);i++)
{
    one_4nc(i) = 1;
}

    //Liquid phase fugacity calculation
    //am and bm calculation
    phase = 1; //1 for liquid, 2 for vapor
    bl = b_mixing_rules_function(nc, b, x, MR);
    al = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    phase = 2;
    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    if(counter == 0 || counter == max_num_iter)
    {
        Vlinit = bl/0.99; //iota == 0.99
        Vvinit = bv+(R*T/P); //iota = bv/(bv+(R*T/P), Vvinit = bv/iota

        if(iter_choice==1)
        {
        Vl = Vlinit;
        Vv = Vvinit;
        }
    }


    if(EdE==3 || EdE==6)
    {
    deltaV = 0;
    phase = 1;
    Xl = volume_function(nc, EdE, phase, x, Xl, EdE_parameters, bl, al, R, T, P, tolV, tolZl, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto,
                        beta_auto, Dij);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto,
                        beta_auto, Dij);
    }
    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);

    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vl1, Vl2);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
Kx = (x.asDiagonal())*K;

    for(i=0; i<nc; i++)
    {
         one(i) = 1;
    }

sumKx = one.transpose()*Kx;
double sumKxold;
sumKxold = sumKx;
errorSUMKx = tolSUMKx + 1;


double counter2;
counter2 = 0;

while(errorSUMKx>tolSUMKx || counter2<=1)
    {
    y = Kx.array()/sumKx;

    initialSUMKx = sumKx;

    //Vapor phase fugacity calculation
    //am and bm calculation
    if(process==2)
    {
    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    phase = 2; //1 for liquid, 2 for vapor

    if(EdE==3 || EdE==6)
    {
    phase = 2;
    deltaV = 0;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto,
                        beta_auto, Dij);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);
    }

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


    K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
    Kx = (x.asDiagonal())*K;
    sumKxnew = one.transpose()*Kx;

    finalSUMKx = sumKxnew;
    errorSUMKx = fabs(finalSUMKx - initialSUMKx);

    errorSUMKx = errorSUMKx/finalSUMKx;

    sumKx = sumKxnew;

 if(counter2==200)
 {
   sumKx = sumKxold;
   errorSUMKx = 0.00000000000001;
 }
 counter2++;

    }

double errorKxnew;
Ey = sumKx-1;
errorKx = fabs(Ey);

errorKx = errorKx/sumKx;

y = Kx.array()/sumKx;


switch(process)
{
case 1: //Isothermic
P = P*sumKx;
break;

case 2: //Isobaric

    Told = T;

    T = 0.1*T/sumKx+0.9*T; //AQUI DIVIDE

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();
break;
}

if(isnan(errorKx)==1 && process==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
}

counter++;

double trivial, V_check;

switch(process)
{
case 1: //Isothermic
    if(isnan(P)==1 || isinf(P)==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "P = NAN OR INF" << endl;
    cin.get();
}
break;

case 2: //Isobaric
    if(isnan(T)==1 || isinf(T)==1)
{
    T = (1+0.1*k)*T;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "T = NAN OR INF" << endl;

    counter = 500;
    //cin.get();
}
break;
}

if(counter==max_num_iter)
    {
    errorKx=0.00000000000001;
    }
}

if(counter!=max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    Pinit = P;
    break;

    case 2: //Isobaric
    Tinit = T;
    break;
    }
    Vlinit = Vl;
    Vvinit = Vv;
    yinit = y;
}

if(counter==max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
        Pinit = Psat.transpose()*x;
        break;

    case 2: //Isobaric
        Tinit = Tsat.transpose()*x;
        break;
    }

    yinit = ((Psat*x.transpose()).diagonal()).array()/Pinit;
}

//------------------------------------------
//P = P*100; //Converting from bar para kPa
//Converting directly on output
y = Kx;
gama = Psat.asDiagonal().inverse()*(x.asDiagonal().inverse()*y);
gama = P*gama.array();
ln_gama = gama.array().log();
G_ex = x.transpose()*ln_gama;
G_ex = G_ex*R*T;
cout << "--------------------------------" << endl;
cout << "Zl = " << Zl << endl;
cout << "Zv = " << Zv << endl;
cout << "Vl = " << Vl << endl;
cout << "Vv = " << Vv << endl;
cout << "x1 = " << x(0) << endl;
cout << "y1 = " << y(0) << endl;
cout << "P(bar) = " << P << endl;
cout << "T(K) = " << T << endl;
cout << "errorKx = " << errorKx << endl;
cout <<"--------- counter = " << counter << " ---------" << endl;
    if(process==1)
    {
    output << x(0) << ";" << y(0) << ";" << T << ";" << P*100 << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex << endl;
    }

    if(process==2)
    {
    output << x(0) << ";" << y(0) << ";" << P*100 << ";" << T << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex << endl;
    }


if(iter_choice==1)
 {
 cout << "End of calculation \n \n";
 counter = 0;
 }


Pold = P;



}

}

}

//CÁLCULO PARA MISTURA BINÁRIA
//==============================================================================================================
if(mixture==2)
{

for (x(0)=0.0001 ; x(0)<=1.000 ; x(0)=x(0)+0.005)
{
 x(1) = 1-x(0);

 if(iter_choice==1)
 {
 cout << "Define x for component 1: ";
 cin >> x(0);
 x(1) = 1-x(0);

 counter = 0;

    if(counter==0)
    {
        switch(process)
        {
            case 1: //Isothermic
            P = Psat.transpose()*x;
            Pinit = P;
            break;

            case 2: //Isobaric
            T = Tsat.transpose()*x;
            Tinit = T;
            CT = C.array()+T;
            logPsat = A - (CT.asDiagonal().inverse()*B);
            ln10.fill(log(10));
            lnPsat = (logPsat*ln10.transpose()).diagonal();
            Psat = lnPsat.array().exp();
            cout << "Tinit" << Tinit << endl;
            break;
            //cin.get();
        }

    y = ((Psat*x.transpose()).diagonal()).array()/P;
    yinit = y;
    }

 }


 if(nc>2)
 {
     iter_choice=1;

     int q;
     q = 0;
     while(q<nc)
     {
     cout << "Define x for component " << q+1 << " : ";
     cin >> x(q);
     q++;
     }
 }

switch(process)
{
case 1: //Isothermic
P = Pinit;
break;

case 2: //Isobaric
T = Tinit;
//ATUALIZAR PSAT??????????????????????????????????????????????????????????????????
CT = C.array()+T;
logPsat = A - (CT.asDiagonal().inverse()*B);
ln10.fill(log(10));
lnPsat = (logPsat*ln10.transpose()).diagonal();
Psat = lnPsat.array().exp();

Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
break;
}

if(counter == max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    P = Psat.transpose()*x;
    break;

    case 2: //Isobaric
    T = Tsat.transpose()*x;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    break;
    }

y = ((Psat*x.transpose()).diagonal()).array()/P;
}

if(counter != max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    P = Pinit;

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }
    break;

    case 2: //Isobaric
    T = Tinit;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();

        if(isnan(y(0))==1 || isinf(y(0))==1)
        {
        y = ((Psat*x.transpose()).diagonal()).array()/P;
        }

        else
        {
        y = yinit;
        }
    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    break;
    }

Vlinit = Vl;
Vvinit = Vv;
}

counter = 0;
int k;
k=1;
errorKx = tolKx + 1;
tol_u = 0.00001;

cout << "x = " << x(0) << endl;

while(errorKx>tolKx)
{
Tr = T*Tc.asDiagonal().inverse().diagonal();
alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);

double Bcpa;
MatrixXd pre_F(nc,4);
    VectorXd one_4(4), one_4nc(4*nc);
    Bcpa = x.transpose()*b;
one_4 <<    1,
            1,
            1,
            1;

for(i=0;i<(4*nc);i++)
{
    one_4nc(i) = 1;
}

    //Liquid phase fugacity calculation
    //am and bm calculation
    phase = 1; //1 for liquid, 2 for vapor
    bl = b_mixing_rules_function(nc, b, x, MR);
    al = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    phase = 2;
    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    if(counter == 0 || counter == max_num_iter)
    {
        Vlinit = bl/0.99; //iota == 0.99
        Vvinit = bv+(R*T/P); //iota = bv/(bv+(R*T/P), Vvinit = bv/iota

        if(iter_choice==1)
        {
        Vl = Vlinit;
        Vv = Vvinit;
        }
    }


    if(EdE==3 || EdE==6)
    {
    deltaV = 0;
    phase = 1;
    Xl = volume_function(nc, EdE, phase, x, Xl, EdE_parameters, bl, al, R, T, P, tolV, tolZl, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto,
                        beta_auto, Dij);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto,
                        beta_auto, Dij);
    }
    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);

    if(EdE==5 || EdE==6)
    {
    phase = 1;
    V_renormalized(phase,x(0),P,bl,R,T,d2P,d2u,x_rv,rho_rv,Pmat,umat,&Vl);
    //V_renormalized(phase,x(0),P,b(0),R,T,d2P1,d2u1,x_rv,rho_rv,P1mat,u1mat,&Vl1);
    //V_renormalized(phase,x(0),P,b(1),R,T,d2P2,d2u2,x_rv,rho_rv,P2mat,u2mat,&Vl2);

    phase = 2;
    V_renormalized(phase,y(0),P,bv,R,T,d2P,d2u,x_rv,rho_rv,Pmat,umat,&Vv);
    //V_renormalized(phase,y(0),P,b(0),R,T,d2P1,d2u1,x_rv,rho_rv,P1mat,u1mat,&Vv1);
    //V_renormalized(phase,y(0),P,b(1),R,T,d2P2,d2u2,x_rv,rho_rv,P2mat,u2mat,&Vv2);
    //cout << "Vl: " << Vl << " " << Vl1 << " " << Vl2 << " " << x(0)*Vl1+x(1)*Vl2 << endl;
    //cout << "Vv: " << Vv << " " << Vv1 << " " << Vv2 << " " << x(0)*Vv1+x(1)*Vv2 << endl;
    }

    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vl1, Vl2);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
Kx = (x.asDiagonal())*K;
cout << "phi:   am:    bm:   P:   x0:   x1:   " << endl;
cout << phi_liquid_phase << " " << al << " " << bl << " " << P << " " << x(0) << " " << x(1) << endl;
cout << phi_vapor_phase << " " << av << " " << bv << " " << P << " " << y(0) << " " << y(1) << endl;
cin >> stop;

    for(i=0; i<nc; i++)
    {
         one(i) = 1;
    }

sumKx = one.transpose()*Kx;
double sumKxold;
sumKxold = sumKx;
errorSUMKx = tolSUMKx + 1;

double counter2;
counter2 = 0;

while(errorSUMKx>tolSUMKx || counter2<=1)
    {
    y = Kx.array()/sumKx;

    initialSUMKx = sumKx;

    //Vapor phase fugacity calculation
    //am and bm calculation

    if(process==2)
    {
    Tr = T*Tc.asDiagonal().inverse().diagonal();
    alfa = alfa_function(EdE, nc, Tr, omega, a0, c1);
    a = a_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, EdE, a0);
    b = b_function(nc, R, EdE_parameters, omega, Tc, Pc, alfa, bCPA, EdE);
    }

    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    phase = 2; //1 for liquid, 2 for vapor

    if(EdE==3 || EdE==6)
    {
    phase = 2;
    deltaV = 0;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto,
                        beta_auto, Dij);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);
    }

    if(EdE==5 || EdE==6)
    {
    phase = 2;
    V_renormalized(phase,y(0),P,bv,R,T,d2P,d2u,x_rv,rho_rv,Pmat,umat,&Vv);
    //V_renormalized(phase,y(0),P,b(0),R,T,d2P1,d2u1,x_rv,rho_rv,P1mat,u1mat,&Vv1);
    //V_renormalized(phase,y(0),P,b(1),R,T,d2P2,d2u2,x_rv,rho_rv,P2mat,u2mat,&Vv2);
    }

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1, d2P, d2u1, d2u2,
                                         x_rv, rho_rv, Pmat, umat, u1mat, u2mat, Vv1, Vv2);


    K = (phi_vapor_phase.asDiagonal().inverse())*phi_liquid_phase;
    Kx = (x.asDiagonal())*K;
    sumKxnew = one.transpose()*Kx;

    finalSUMKx = sumKxnew;
    errorSUMKx = fabs(finalSUMKx - initialSUMKx);

    errorSUMKx = errorSUMKx/finalSUMKx;

    sumKx = sumKxnew;

 if(counter2==100)
 {
   sumKx = sumKxold;
   errorSUMKx = 0.0000000001;
 }
 counter2++;

    }

double errorKxnew;
Ey = sumKx-1;
errorKx = fabs(Ey);

errorKx = errorKx/sumKx;

y = Kx.array()/sumKx;
//cout << "errorKx: " << errorKx << " sumKx: " << sumKx << " K: " << K << endl;


switch(process)
{
case 1: //Isothermic
P = P*sumKx;
break;

case 2: //Isobaric

    Told = T;

    T = 0.1*T/sumKx+0.9*T;

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();
break;
}

if(isnan(errorKx)==1 && process==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
}

counter++;

double trivial, V_check;

switch(process)
{
case 1: //Isothermic
    if(isnan(P)==1 || isinf(P)==1)
{
    P = (1+0.1*k)*Pinit;
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "P = NAN OR INF" << endl;
    cin.get();
}
break;

case 2: //Isobaric
    if(isnan(T)==1 || isinf(T)==1)
{
    T = (1+0.1*k)*T;
    CT = C.array()+T;
    logPsat = A - (CT.asDiagonal().inverse()*B);
    ln10.fill(log(10));
    lnPsat = (logPsat*ln10.transpose()).diagonal();
    Psat = lnPsat.array().exp();
    y = ((Psat*x.transpose()).diagonal()).array()/P;
    errorKx = 1;
    k++;
    cout << "T = NAN OR INF" << endl;

    counter = 500;

    //cin.get();
}
break;
}

if(counter==max_num_iter)
    {
    errorKx=0.00000000000001;
    }
}



if(counter!=max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
    Pinit = P;
    break;

    case 2: //Isobaric
    Tinit = T;
    break;
    }
    Vlinit = Vl;
    Vvinit = Vv;
    yinit = y;
}

if(counter==max_num_iter)
{
    switch(process)
    {
    case 1: //Isothermic
        Pinit = Psat.transpose()*x;
        break;

    case 2: //Isobaric
        Tinit = Tsat.transpose()*x;
        break;
    }

    yinit = ((Psat*x.transpose()).diagonal()).array()/Pinit;
}

//P = P*100; //Converting from bar para kPa
//Converting directly on output
y = Kx;
gama = Psat.asDiagonal().inverse()*(x.asDiagonal().inverse()*y);
gama = P*gama.array();
ln_gama = gama.array().log();
G_ex = x.transpose()*ln_gama;
G_ex = G_ex*R*T;
cout << "--------------------------------" << endl;
cout << "Zl = " << Zl << endl;
cout << "Zv = " << Zv << endl;
cout << "Vl = " << Vl << endl;
cout << "Vv = " << Vv << endl;
cout << "x1 = " << x(0) << endl;
cout << "y1 = " << y(0) << endl;
cout << "P(bar) = " << P << endl;
cout << "T(K) = " << T << endl;
cout << "errorKx = " << errorKx << endl;
cout <<"--------- counter = " << counter << " ---------" << endl;

    if(process==1)
    {
    output << x(0) << ";" << y(0) << ";" << T << ";" << P*1000 << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex << endl;
    }

    if(process==2)
    {
    output << x(0) << ";" << y(0) << ";" << P*1000 << ";" << T << ";" << Vl << ";" << Vv << ";"
           << sumKx << ";" << counter << ";" << u_liquid1 << ";" << u_vapor1 << ";"
           << X1l << ";" << X1v << ";" << Zl << ";" << Zv << ";" << phi_liquid_phase(0)
           << ";" << phi_liquid_phase(1) << ";" << phi_vapor_phase(0) << ";" << phi_vapor_phase(1)
           << ";" << Vl_obj << ";" << Vv_obj << ";" << dP_dVl << ";" << dP_dVv << ";" << G_ex << endl;
    }


if(iter_choice==1)
 {
 cout << "End of calculation \n \n";
 counter = 0;
 }

}

}

}

}
