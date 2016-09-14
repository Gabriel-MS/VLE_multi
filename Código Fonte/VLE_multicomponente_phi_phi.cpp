/*
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                ||PROGRAMA PARA CÁLCULO DE EQUILÍBRIO VLE EM MISTURAS MULTICOMPONENTES||
                ||AUTOR: GABRIEL MORAES SILVA                                         ||
                ||LINGUAGEM: C++                                                      ||
                ||BIBLIOTECAS: EIGEN                                                  ||
                ||ANO: 2016                                                            ||
                ||VERSÃO 1.0                                                          ||
                ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/


#include <cmath>
#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>2

typedef std::numeric_limits< double > dbl;

#include "mixingrules.h"
#include "EdE.h"
#include "Gibbs.h"
#include "Association.h"
#include "Renormalization.h"
#include "interpolation_util.h"
#include "numerical.h"

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

//Arquivo de saída dos dados--------------------------------------------------------------------------
ofstream output("../Planilhas de análise/Output.csv");
output << "Dados da simulação \n -------------------------------------------------------------------------------------" << endl;


//Apresentação e versão do programa
cout << "    |===========================||" << endl;
cout << "    |Autor: Gabriel Moraes Silva||" << endl;
cout << "    |Ano: 2016                  ||" << endl;
cout << "    |V. 1.0                     ||" << endl;
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
double E_a[nc], beta_a[nc], a0_a[nc], bCPA_a[nc], c1_a[nc];
double Tc_v[nc], Pc_v [nc], omega_v[nc];
double P, T, R, Pnew, sumKx, al, bl, av, bv, Ey;
double tolZv, tolZl, tolSUMKx, tolKx, initialSUMKx, sumKxnew, finalSUMKx, tolX, tolV;
double errorZv, errorZl, errorSUMKx, errorKx, k12, V, dP_dV, rho_l, X1, Vl, Vv, Vt, deltaV;
double Pinit, Vlinit, Vvinit, Zl, Zv, X1l, X1v;
VectorXd ln_Ki(nc), pre_P1(nc), Tr_1(nc), pre_P1_exp(nc), pre_ln_Ki;
double Vl_obj, Vv_obj, Ql, Qv, dP_dVl, dP_dVv;
double log10P, Tb, Tinit, Told;
double G_ex;
VectorXd Tsat(nc), Alog10P(nc), gama(nc), ln_gama(nc);
double init_T, final_T, step, BETCR, init_P, final_P, Pold;
int max_num_iter, counter, stop, Renormalization;
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
cout << "\nChoose the EoS: \n 1.Soave-Redlich-Kwong \n 2.Peng-Robinson \n 3.CPA-SRK " << endl;
cin >> EdE;
output << "Equation of state = " << EdE << endl;

cout << "\nConsider Renormalization? \n 1.Yes \n 2.No " << endl;
cin >> Renormalization;
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
}

cout << "Tc = " << Tc << endl;
cout << "Tc_v = " << Tc_v << endl;
cout << "Pc_v = " << Pc_v << endl;
cout << "omega_v = " << omega_v << endl;

//Reading C++ vectors into Eigen type vectors
for(n=0; n<nc; n++)
{

if(EdE==3)
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
}

cout << "Tc = " << Tc << endl;
cout << "Tc_virtual = " << Tc_virtual << endl;
cout << "Pc_virtual = " << Pc_virtual << endl;
cout << "omega_virtual = " << omega_virtual << endl;

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
}

}


//Ideal gas constant
R = 0.08314462; // L.bar/K/mol

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
cout << "Psat = " << Psat << endl;
break;

case 2: //Isobaric
log10P = log10(P);
Alog10P = A.array()-log10P;
Tsat = (Alog10P.asDiagonal().inverse()*B)-C;
cout << "Tsat = " << Tsat << endl;
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

if(EdE==3)
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
cout << "y initial guess = \n" << y << endl;
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

    ofstream Renorm("../Planilhas de análise/Renormalization.csv");
    //Renorm << "rho" << ";" << "f" << ";" << "f0" << endl;
    Renorm << "Density" << ";" << "f" << ";" << "f0" << ";" << "Chem. Pot." << ";" << "P" << ";" << "T" << endl;

    int i, j, n;
    long double kB, L, L3, fi, K, rho, rho_plus, rho_minus, fl_plus, fl, fl_minus, fs_plus, fs, fs_minus;
    long double Gl, Gs, OMEGA, delta_f, f, f0, fl_old_plus, fl_old_minus, fs_old_plus, fs_old_minus, fl_old, fs_old;
    long double OMEGAs, OMEGAl, f_old, alfa_r, am, rho_max, bm, tolZ, rho2, var, f0_plus, f0_minus;
    long double width, suml, sums, m, fl_plus_old, fl_minus_old, fs_plus_old, fs_minus_old, f_original;
    long double Gl0, Gs0, Gln, Gsn, eGl0, eGs0, eGln, eGsn, phi_r, P_test_old, P_average_0;
    long double pmax_cond, P_max, P_min, P_average, P_test, test, P_l, u_l, P_v, u_v, pmin_cond, rho_v;
    std::vector<double> rho_vec(1000), f_vec(1000), u_vec(1000), P_vec(1000), f0_vec(1000), dP_dV(1000);

    double Q_func;

    VectorXd L3_v(nc), L_v(nc), fi_v(nc), lnXX2(4*nc), f_assoc(nc), one_4c(4*nc, 2);

    VectorXd rhov(500), fv(500), x_(500), fv_(500);

    n = 100;

    VectorXd fl_old_p(n), fl_oldv(n), fl_old_m(n), fs_old_p(n), fs_oldv(n), fs_old_m(n);

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
    L_v << 7.40e10, 7.50e10;
    fi_v << 8.39, 8.98;

    bm = b_mixing_rules_function(nc, b, x, MR);
    am = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);
    alfa_r = am;


    kB = 1.38064852e-25; //Boltzmann constant L.bar/K

    //L3_v = L_v.array().pow(3);
    //L3 = x.transpose()*L3_v;
    //L = pow(L3,(1/3)); //Cut-off Length

    if(EdE==1)
    {
    if(cp[0]==12)
    {
        L = 7.80e-9;
    }

    else
    {
        cout << "L = " << endl;
        cin >> L;
    }
    }


    if(EdE==3)
    {
    if(cp[0]==47)
    {
        L = 5.51e-09;
        phi_r = 8.39;
    }

    else
    {
        cout << "L = " << endl;
        cin >> L;
        cout << "phi_r = " << endl;
        cin >> phi_r;
    }
    }
    //fi = x.transpose()*fi_v; //Second crossover parameter phi

    rho_max = 0.99999/bm;


    int t, k, w;
    t = 0;
    k = 0;

    cout << "\nDefine final Temperature: ";
    cin >> final_T;

    cout << "\nDefine steps: ";
    cin >> step;

    init_T = T;
    Told = T;

//while(T<(Tc[0]+10))
for(T=init_T; T<=final_T; T=T+step)
{

    E_row = ((E_row.array().log())*Told/T).exp();
    E_col = ((E_col.array().log())*Told/T).exp();
    E_auto = ((E_auto.array().log())*Told/T).exp();

    Told = T;

    rho = 0.0001;
    w = 0;

while(rho<=rho_max)
{
    V = 1/rho;

    if(EdE==3)
    {
    X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, V, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);
    }

    rho2 = min(rho, rho_max-rho);

    f0 = helmholtz_repulsive(EdE, R, T, rho, am, bm, X, x);
    f_original = f0;

    f0 = f0 + 0.5*am*rho*rho;
    f_old = f0;

    for(i=1;i<11;i=i+1)
    {
        K = kB*T/((pow(2,3*i))*pow(L,3));

            if(i==1)
            {
            fl_old = f_old;
            fs_old = f_old;
            }

        width = (rho2-0)/n;
        suml = 0;
        sums = 0;
        m = n-1;


        for(j=0;j<n;j++)
        {

            var = 0+j*width;
            rho_plus = rho+var;
            rho_minus = rho-var;
/*
            if(rho_minus <= 0)
            {
                rho_minus = 0;
            }

            if(rho_plus >= rho2)
            {
                rho_plus = rho2;
            }
*/
            if(i==1)
            {
            fl_old_p(j) = 0.5*am*(rho_plus)*(rho_plus)+helmholtz_repulsive(EdE, R, T, rho_plus, am, bm, X, x);
            fl_oldv(j) = 0.5*am*(rho)*(rho)+helmholtz_repulsive(EdE, R, T, rho, am, bm, X, x);
            fl_old_m(j) = 0.5*am*(rho_minus)*(rho_minus)+helmholtz_repulsive(EdE, R, T, rho_minus, am, bm, X, x);

            fs_old_p(j) = 0.5*am*(rho_plus)*(rho_plus)+helmholtz_repulsive(EdE, R, T, rho_plus, am, bm, X, x);
            fs_oldv(j) = 0.5*am*(rho)*(rho)+helmholtz_repulsive(EdE, R, T, rho, am, bm, X, x);
            fs_old_m(j) = 0.5*am*(rho_minus)*(rho_minus)+helmholtz_repulsive(EdE, R, T, rho_minus, am, bm, X, x);
            }

/*
            if(i==1)
            {
            fl_old_p(j) = helmholtz_repulsive(EdE, R, T, rho_plus, am, bm, X, x);
            fl_oldv(j) = helmholtz_repulsive(EdE, R, T, rho, am, bm, X, x);
            fl_old_m(j) = helmholtz_repulsive(EdE, R, T, rho_minus, am, bm, X, x);

            fs_old_p(j) = helmholtz_repulsive(EdE, R, T, rho_plus, am, bm, X, x);
            fs_oldv(j) = helmholtz_repulsive(EdE, R, T, rho, am, bm, X, x);
            fs_old_m(j) = helmholtz_repulsive(EdE, R, T, rho_minus, am, bm, X, x);
            }
*/

            fl_plus = helmholtz_recursion_long(EdE, fl_old_p(j), rho_plus, am);
            fl = helmholtz_recursion_long(EdE, fl_oldv(j), rho, am);
            fl_minus = helmholtz_recursion_long(EdE, fl_old_m(j), rho_minus, am);

            fs_plus = helmholtz_recursion_short(EdE, fs_old_p(j), rho_plus, am, i, L, phi_r);
            fs = helmholtz_recursion_short(EdE, fs_oldv(j), rho, am, i, L, phi_r);
            fs_minus = helmholtz_recursion_short(EdE, fs_old_m(j), rho_minus, am, i, L, phi_r);

            //if(i!=1)
            //{
            fl_old_p(j) = fl_plus;
            fl_oldv(j) = fl;
            fl_old_m(j) = fl_minus;

            fs_old_p(j) = fs_plus;
            fs_oldv(j) = fs;
            fs_old_m(j) = fs_minus;
            //}
/*
            fl_plus = helmholtz_recursion_long(EdE, fl_old, rho+var, am);
            fl = helmholtz_recursion_long(EdE, fl_old, rho, am);
            fl_minus = helmholtz_recursion_long(EdE, fl_old, rho-var, am);

            fs_plus = helmholtz_recursion_short(EdE, fs_old, rho+var, am, i, L, phi_r);
            fs = helmholtz_recursion_short(EdE, fs_old, rho, am, i, L, phi_r);
            fs_minus = helmholtz_recursion_short(EdE, fs_old, rho-var, am, i, L, phi_r);
*/
            Gl = (fl_plus-2*fl+fl_minus)/2;
            Gs = (fs_plus-2*fs+fs_minus)/2;

            if(j==0 || j==(n-1))
            {
            suml = suml + 0.5*exp(-Gl/K);
            sums = sums + 0.5*exp(-Gs/K);
            }

            else
            {
            suml = suml + exp(-Gl/K);
            sums = sums + exp(-Gs/K);
            }

            //cout << "j = " << j << " // fl_plus = " << fl_plus << " // fl = " << fl << " // fs_plus = "
            //<< fs_plus << " // Gl = " << Gl << " // Gs = " << Gs << endl;

            //cout << "j = " << j << "// rho_minus = " << rho_minus << "// rho_max = " << rho_plus << endl;

        }
            //fl_old = fl;
            //fs_old = fs;

        OMEGAl = width*suml;
        OMEGAs = width*sums;

        delta_f = -K*log(OMEGAs/OMEGAl);

        if(rho>=(rho_max/2))
        {
            delta_f = 0;
        }

        f = f_old + delta_f;

        f_old = f;
        //fl_old = f_old;
        //fs_old = f_old;

        //cout << "i = " << i << " // f = " << f << " // delta_f = " << delta_f <<
        //" // Is = " << OMEGAs << " // Il = " << OMEGAl << endl;
    }
    f = f - 0.5*am*rho*rho;

    rho_vec[w] = rho;
    f_vec[w] = f;
    f0_vec[w] = f0;

    //cout << "rho = " << rho << "  //  f = " << f << endl;
    //Renorm << std::fixed << std::setprecision(15) << rho << ";" << f << ";" << f_original << ";" << T << endl;

    rho = rho+rho_max/1000;

    w++;
}

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec);

for(i=0; i<1000; i++)
{
    P_vec[i] = -f_vec[i] + rho_vec[i]*u_vec[i];
    //cout << P_vec[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec[i] << ";" << f_vec[i] << ";"
    //       << f0_vec[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << T << endl;
}


//Search Max and Min Pressures from isotherm
//

dP_dV[0] = P_vec[0]/rho_vec[0]; //Wrong, fix

  for(i=1; i<1000; i++)
  {
    dP_dV[i] = (P_vec[i]-P_vec[i-1])/(rho_vec[i]-rho_vec[i-1]);
  }

//Maximum pressure at dP/dV=0
i=0;
do
{
 pmax_cond = dP_dV[i];
 i++;
}while(pmax_cond>0);
P_max = P_vec[i-2];
cout << "dP/dV at max = " << dP_dV[i-2] << endl;

//Minimum pressure at dP/dV=0
j=999;
do
{
 pmin_cond = dP_dV[j];
 j--;
//cout << "dP/dV = " << dP_dV[j] << endl;
}while(pmin_cond>0);
P_min = P_vec[j+2];
cout << "dP/dV at max = " << dP_dV[i+2] << endl;

cout << "max pressure at = " << P_max << endl;
cout << "min pressure at = " << P_min << endl;
//
//------------------------------------------

errorKx = 1e10; //Forces to enter loop (could use "do" command, but i don't want)
j=0;
counter = 0;

int tl, tv, c;
c = 1;

//THIS SHOULD BE INSIDE ERRORKX ITERATION
double tol_P = 1e-3;

if(counter==0)
{
P_average = (P_max+P_min)/2;
P_average_0 = P_average;
cout << "P_average = " << P_average << endl;
}

cout << "Pressure_average = " << P_average << endl;


//Liquid phase
i=0;
P_test_old = 1e10;
do
{
    P_test = fabs(P_average-P_vec[i]);
    i++;
    //cout << "P_test_l = " << P_test << endl;

    if(P_test>P_test_old) P_test = 1e-4; // Guarantees getting the nearest value

    P_test_old = P_test;
}while(P_test > tol_P);
P_v = P_vec[i-1];
u_v = u_vec[i-1];
rho_v = rho_vec[i-1];
Vv = 1/rho_v;
tv = i-1;

//Vapor phase
j=999;
P_test_old = 1e10;
do
{
    P_test = fabs(P_average-P_vec[j]);
    j--;
    //cout << "P_test_v = " << P_test << endl;

    if(P_test>P_test_old) P_test = 1e-4; // Guarantees getting the nearest value

    P_test_old = P_test;
}while(P_test > tol_P);
P_l = P_vec[j+1];
u_l = u_vec[j+1];
rho_l = rho_vec[j+1];
Vl = 1/rho_l;
tl = j+1;
//END ERRORKX ITERATION PART

cout << "P_l = " << P_l << endl;
cout << "P_v = " << P_v << endl;
cout << "rho_l = " << rho_l << endl;
cout << "rho_v = " << rho_v << endl;


/*
while(errorKx > tolKx)
{
//Choose medium value of Pressure, get nearest P in P_vec and select rho_L and rho_V accordingly with chemical potentials
//
double tol_P = 1e-3;

if(counter==0)
{
P_average = (P_max+P_min)/2;
P_average_0 = P_average;
cout << "P_average = " << P_average << endl;
}

cout << "Pressure_average = " << P_average << endl;

//Liquid phase
i=0;
P_test_old = 1e10;
do
{
    P_test = fabs(P_average-P_vec[i]);
    i++;
    //cout << "P_test_l = " << P_test << endl;

    if(P_test>P_test_old) P_test = 1e-4; // Guarantees getting the nearest value

    P_test_old = P_test;
}while(P_test > tol_P);
P_l = P_vec[i-1];
u_l = u_vec[i-1];
rho_l = rho_vec[i-1];
Vl = 1/rho_l;
tl = i-1;

//Vapor phase
j=999;
P_test_old = 1e10;
do
{
    P_test = fabs(P_average-P_vec[j]);
    j--;
    //cout << "P_test_v = " << P_test << endl;

    if(P_test>P_test_old) P_test = 1e-4; // Guarantees getting the nearest value

    P_test_old = P_test;
}while(P_test > tol_P);
P_v = P_vec[j+1];
u_v = u_vec[j+1];
rho_v = rho_vec[j+1];
Vv = 1/rho_v;
tv = j+1;

cout << "P_l = " << P_l << endl;
cout << "P_v = " << P_v << endl;
cout << "rho_l = " << rho_l << endl;
cout << "rho_v = " << rho_v << endl;
//
//------------------------------------------

//Calculate phi_liquid_phase and vapor_phase
//
    phase = 1; //1 for liquid, 2 for vapor
    bl = b_mixing_rules_function(nc, b, x, MR);
    al = a_mixing_rules_function(nc, MR, a, x, k12, bl, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    phase = 2;
    bv = b_mixing_rules_function(nc, b, y, MR);
    av = a_mixing_rules_function(nc, MR, a, y, k12, bv, b, T, q, r, Aij, R, alfa_NRTL, EdE, G_ex_model);

    if(EdE==3)
    {
    deltaV = 0;

    phase = 1;
    Xl = fraction_nbs(nc, combining_rule, phase, R, T, P_l, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, Vl, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);
    phase = 2;
    Xv = fraction_nbs(nc, combining_rule, phase, R, T, P_v, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, y, EdE, EdE_parameters, b, tolZ, Vv, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);
    }

    phi_liquid_phase[0] = (u_l-R*T*log(P_l))/R/T-P_l/(R*T*rho_l);
    phi_vapor_phase[0] = (u_v-R*T*log(P_v))/R/T-P_v/(R*T*rho_v);
//
//------------------------------------------

errorKx = fabs(phi_liquid_phase[0]/phi_vapor_phase[0]-1);
cout << "errorKx = " << errorKx << endl;

//Next value of pressure, or get out the loop
//

if(c%2==0)
{
P_average = P_average_0*(1+c/100);
}

else
{
P_average = P_average_0-(1-c/100);
}
counter++;
j++;
c++;
cout << "c = " << c << endl;
//
//------------------------------------------
}
*/

cout << "T = " << T << endl;

k++;

Renorm << std::fixed << std::setprecision(15) << T << ";" << rho_l << ";"
           << rho_v << ";" << P_l << ";" << P_v << ";" << endl;

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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto, beta_auto);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);

    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);

    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);


    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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


    if(EdE==3)
    {
    deltaV = 0;
    phase = 1;
    Xl = volume_function(nc, EdE, phase, x, Xl, EdE_parameters, bl, al, R, T, P, tolV, tolZl, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto, beta_auto);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);
    }
    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);

    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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

    if(EdE==3)
    {
    phase = 2;
    deltaV = 0;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);
    }

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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


    if(EdE==3)
    {
    deltaV = 0;
    phase = 1;
    Xl = volume_function(nc, EdE, phase, x, Xl, EdE_parameters, bl, al, R, T, P, tolV, tolZl, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR, E_auto, beta_auto);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);
    }
    X1l = Xl(0)*Xl(1)*Xl(2)*Xl(3);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);
    phase = 1;
    phi_liquid_phase = fugacity_function(nc, phase, al, bl, a, b, R, T, P, tolZl, EdE_parameters, MR, q_prime, r, Aij, x, q, EdE,
                                         alfa_NRTL, G_ex_model, k12, Xl, tolV, Vl, n_v, Vl, &Zl, &u_liquid1);

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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

    if(EdE==3)
    {
    phase = 2;
    deltaV = 0;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR, E_auto, beta_auto);
    X1v = Xv(0)*Xv(1)*Xv(2)*Xv(3);
    }

    phase = 2;
    phi_vapor_phase = fugacity_function(nc, phase, av, bv, a, b, R, T, P, tolZv, EdE_parameters, MR, q_prime, r, Aij, y, q, EdE,
                                        alfa_NRTL, G_ex_model, k12, Xv, tolV, Vv, n_v, Vv, &Zv, &u_vapor1);


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

}

}

}

}
