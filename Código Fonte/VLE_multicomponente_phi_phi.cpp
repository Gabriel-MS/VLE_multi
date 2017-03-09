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
int max_num_iter, counter, stop, Renormalization, sr_type, r_type, Iteration;
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

cout << "\nConsider Renormalization? \n 1.Yes \n 2.No " << endl;
cin >> Renormalization;

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
}

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
R = 0.000008314; //MPa.m³/mol/K
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

    ofstream Renorm;
    ofstream Not_splined;
    ofstream Envelope;
    ofstream dfnout;
    //ofstream before_renorm;
    //ofstream after_renorm;

    Renorm.open("../Planilhas de análise/Renormalization.csv");
    Not_splined.open("../Planilhas de análise/Renormalization_not_splined.csv");
    dfnout.open("dfn_msa_out.csv");
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

    int i, j, n, k, w, t;
    long double kB, L, L3, fi, K, rho, rho_plus, rho_minus, fl_plus, fl, fl_minus, fs_plus, fs, fs_minus;
    long double Gl, Gs, OMEGA, delta_f, f, f0, fl_old_plus, fl_old_minus, fs_old_plus, fs_old_minus, fl_old, fs_old;
    long double OMEGAs, OMEGAl, f_old, alfa_r, am, rho_max, bm, tolZ, rho2, var, f0_plus, f0_minus;
    long double width, suml, sums, m, fl_plus_old, fl_minus_old, fs_plus_old, fs_minus_old, f_original;
    long double Gl0, Gs0, Gln, Gsn, eGl0, eGs0, eGln, eGsn, phi_r, P_test_old, P_average_0;
    long double pmax_cond, P_max, P_min, P_average, P_test, test, P_l, u_l, P_v, u_v, pmin_cond, rho_v;
    double pi, eps, lambda, sigma, zeta_squared, NA;
    double sig1, sig2, sig12, l12, cnst;

    std::vector<double> rho_vec_out(1000), dP2dV2(1000), dP_dV(1000), P_vec(1000), du_dV(1000);
    std::vector<double> u_vec(1000), u_vec_0(1000), P_vec_0(1000), f_vec_out(1000), f0_vec_out(1000);


    //MatrixXd Area(1000,1000);

    double Q_func, Kn, Ins, Inl, aminl, amins, al, as;
    double Area;
    int flag = 0;

    VectorXd L3_v(nc), L_v(nc), fi_v(nc), lnXX2(4*nc), f_assoc(nc), one_4c(4*nc, 2);

    VectorXd rhov(500), x_(500), fv_(500), X_plus(4*nc), X_minus(4*nc);

    cout << "Renormalization density steps: ";
    cin >> n;
    //n = 5000;


    std::vector<double> rho_vec(n), f_vec(n), u_vec1(n), f0_vec(n), P_vec1(n), Glv2(n), Gsv2(n);
    std::vector<double> flvv(n), fsvv(n);
    VectorXd fl_old_p(n), fl_oldv(n), fl_old_m(n), fs_old_p(n), fs_oldv(n), fs_old_m(n), rho_vector2(n), f_after(n), f_before(n);
    VectorXd flv(n), fsv(n), fv(n), rho_vector(n), delta_fv(n), f_originalv(n), Glv(n), Gsv(n), argl(n), args(n);
    VectorXd P_vec_e(1000), u_vec_e(1000), rho_vec_out_e(1000), f_vec_e(1000);
    std::vector<double> V_vec(1000), A_vec(1000), f_env(1000), rho_env(1000), P_env(1000), u_env(1000);

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
if(EdE != 4)
{
        cout << "\nL (m) = ";
        cin >> L;
        cout << endl;

        cout << "\nphi = ";
        cin >> phi_r;
        cout << endl;
}
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
    cout << "Type of envelope: \n1.Manual \n2.Automatic \n";
    cin >> critical_find;

    /*
    int estimation;
    cout << "Adjust to experimental data estimating L and phi?: \n1.Yes \n2.No \n";
    cin >> estimation;
    */

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


    for(k=0; k<n; k++)
    {
    rho_vec[k] = double(k)/n/bm;
    rho_vector(k) = double(k)/n/bm;
    rho_vec[0] = 1e-6;
    rho_vector(0) = 1e-6;

    //NON BONDED FRACTION DIMENSIONAL**************************
    T = T/bm/R*am;

    if(EdE==3) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_vector(k), deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);

    T = T*bm*R/am;
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    rho_vec[k] = rho_vec[k]*bm;
    rho_vector(k) = rho_vector(k)*bm;
    rho_vec[0] = rho_vec[0]*bm;
    rho_vector(0) = rho_vector(0)*bm;

    fv(k) = helmholtz_repulsive(EdE, R, T, rho_vector(k), am, bm, X, x, sigma, eps, kB);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << endl;

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

        //cout << "Kn = " << Kn << endl;
        //cout << "fv before = " << fv(201) << endl;
        //Kn = kB*T*NA/((pow(2,3*i))*L*b(0));

        //Preparar vetores f_l e f_s
        rho = 1e-6;

        //DIMENSIONLESS!!!************************************************************
        rho = 1e-6*bm;

        if(EdE==3) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);

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


for(w=0; w<1000; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/1000/bm;
    rho_vec_out[0] = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

    if(EdE!=4) T = T*am/bm/R;
    if(EdE==4) T = T*am/R;

for(w=0; w<n; w++)
{
    f_vec[w] = f_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
    f0_vec[w] = f0_vec[w] - rho_vec[w]*R*T*(log(rho_vec[w])-1);
}

f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);

//Add ideal gas contribution before cubic spline

for(w=0; w<n; w++)
{
    f_vec_out[w] = f_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    f0_vec_out[w] = f0_vec_out[w] + rho_vec_out[w]*R*T*(log(rho_vec_out[w])-1);
    u_vec[w] = u_vec[w] + R*T*(log(rho_vec_out[w]));
    u_vec_0[w] = u_vec_0[w] + R*T*(log(rho_vec_out[w]));
}
    if(EdE!=4) T = T/am*bm*R;
    if(EdE==4) T = T/am*R;

for(i=0; i<1000; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b;
    a = (P_vec[505]-P_vec[495])/(rho_vec_out[505]-rho_vec_out[495]);
    b = P_vec[505]-a*rho_vec_out[505];
    P_vec[496] = a*rho_vec_out[496] + b;
    P_vec[497] = a*rho_vec_out[497] + b;
    P_vec[498] = a*rho_vec_out[498] + b;
    P_vec[499] = a*rho_vec_out[499] + b;
    P_vec[500] = a*rho_vec_out[500] + b;
    P_vec[501] = a*rho_vec_out[501] + b;
    P_vec[502] = a*rho_vec_out[502] + b;
    P_vec[503] = a*rho_vec_out[503] + b;
    P_vec[504] = a*rho_vec_out[504] + b;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<1000; i++)
{
// Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
//        << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;

        //DIMENSIONLESS!!!************************************************************
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

if(EdE==4) cout << "T = " << T*am/R << endl;
if(EdE!=4) cout << "T = " << T*am/bm/R << endl;
cout << "=======================================\n" << endl;

k++;

        T = T + step;
        g++;
}

    envelope_tracer(1e-5,env_type);

    break;

 //===============================================================================================================//
 //
 //                                       AUTOMATIC CRITICAL POINT BELOW
 //
 //==============================================================================================================//

    case 2: //Find critical point automatically
    ofstream Envelope("../Planilhas de análise/env.csv");
    ofstream envelope_exponent("../Planilhas de análise/env_exponent.csv");
    ofstream crit_isot("../Planilhas de análise/crit_isotherm.csv");
    Envelope << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << ";" << "delta_P" << ";" << "delta_u" << endl;
    envelope_exponent << "T" << ";" << "rho1" << ";" << "rho2" << ";" << "P" << ";" << "u" << endl;
    double Tcritical, Pcritical, rhocritical1, rhocritical2, rhocritical, ucritical;
    int counter;
    double Tnew;
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

    if(EdE==3) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho_vector(k), deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);

    T = T*bm*R/am;
    //*********************************************************


    //DIMENSIONLESS!!!************************************************************
    rho_vec[k] = rho_vec[k]*bm;
    rho_vector(k) = rho_vector(k)*bm;
    rho_vec[0] = rho_vec[0]*bm;
    rho_vector(0) = rho_vector(0)*bm;

    fv(k) = helmholtz_repulsive(EdE, R, T, rho_vector(k), am, bm, X, x, sigma, eps, kB);
    //cout << "rho / fv = " << rho_vec[k] << " / " << fv(k) << " / "  << X(0) << endl;

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

        if(EdE==3) X = fraction_nbs(nc, combining_rule, phase, R, T, P, tolV, alfa, am, bm, beta_col, beta_row, E_col, E_row,
                     tolX, x, EdE, EdE_parameters, b, tolZ, 1/rho, deltaV, X, 0, a, &Q_func, BETCR, E_auto, beta_auto);

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


for(w=0; w<1000; w++) //FAAAAAAAAAAAAAAAALSOOOOOOOOO
{
    rho_vec_out[w] = double(w)/1000/bm;
    rho_vec_out[0] = 1e-6;

    //DIMENSIONLESS!!!************************************************************
    //rho_vec_out[w] = double(w)/1000/bm*bm;
    //rho_vec_out[0] = 1e-6*bm;
}

f_vec_out = cspline_vec(rho_vec, f_vec, rho_vec_out);
f0_vec_out = cspline_vec(rho_vec, f0_vec, rho_vec_out);

u_vec = cspline_deriv1_vec(rho_vec, f_vec, rho_vec_out);
u_vec_0 = cspline_deriv1_vec(rho_vec, f0_vec, rho_vec_out);


for(i=0; i<1000; i++)
{
    P_vec[i] = -f_vec_out[i] + rho_vec_out[i]*u_vec[i];
    P_vec_0[i] = -f0_vec_out[i] + rho_vec_out[i]*u_vec_0[i];

    //Renorm << std::fixed << std::setprecision(15) << rho_vec_out[i] << ";" << f_vec_out[i] << ";"
    //       << f0_vec_out[i] << ";" << u_vec[i] << ";" << P_vec[i] << ";" << u_vec_0[i] << ";" << P_vec_0[i] << ";" << T << endl;
}

    double a, b;
    a = (P_vec[505]-P_vec[495])/(rho_vec_out[505]-rho_vec_out[495]);
    b = P_vec[505]-a*rho_vec_out[505];
    P_vec[496] = a*rho_vec_out[496] + b;
    P_vec[497] = a*rho_vec_out[497] + b;
    P_vec[498] = a*rho_vec_out[498] + b;
    P_vec[499] = a*rho_vec_out[499] + b;
    P_vec[500] = a*rho_vec_out[500] + b;
    P_vec[501] = a*rho_vec_out[501] + b;
    P_vec[502] = a*rho_vec_out[502] + b;
    P_vec[503] = a*rho_vec_out[503] + b;
    P_vec[504] = a*rho_vec_out[504] + b;
    //cout << "a = " << a << " / " << b << P_vec[500] << endl;

for(i=0; i<1000; i++)
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


for(i=0;i<1000;i++)
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
        dta = dens_newt(rho_vec_out,f_vec_out,P_vec,u_vec,1e-5);

        cout << "\n" << dta[0] << " / " << dta[1] << " / " << dta[2] << " / " << dta[3] << " / " << dta[4] << endl;
        T = T/bm/R*am;
        if(EdE==4) T = T*bm;
        cout << "T = " << T << endl;
        cout << "=======================================\n" << endl;

        k++;
        g++;

        if(flag==6)
        {
        envelope_exponent << T << ";" << dta[0] << ";" << dta[1] << ";" << dta[2] << ";" << dta[4] << endl;
        }

        if(isnan(dta[1])==0)
        {
        Envelope << T << ";" << dta[0] << ";" << dta[1] << ";" << dta[2] << ";" << dta[4] << ";" << fabs(dta[3]-dta[2])
             << ";" << fabs(dta[4]-dta[5]) << endl;
        }

        drho_new = fabs(dta[0]-dta[1]);

        if(counter!=0)
        {
        if(drho_new>drho_old) flag=1;
        }
        counter++;

        if(EdE==3 && flag<5)
        {
        out = T_tracer_CPA(T,drho_new,flag,dta[1]);
        T = out[0];
        flag = out[1];
        }

        if(EdE!=3 && flag<5)
        {
        out = T_tracer(T,drho_new,flag,dta[1]);
        T = out[0];
        flag = out[1];
        }

        drho_old = drho_new;

        Tnew = T;

        if(flag==6)
        {
        cout << "\nflag 6" << endl;
        if(T>=Tcritical) flag=7;
        cout << "T / Tcritical / flag = " << T << " / " << Tcritical << " / " << flag << endl;
            if(T==Tcritical)
            {
            Pcritical = dta[2];
            rhocritical1 = dta[0];
            rhocritical2 = dta[1];
            rhocritical = (rhocritical1+rhocritical2)/2;
            ucritical = dta[4];
            }
        T = T + 0.1;
        Tnew = T;
        }

        if(flag==5)
        {
        cout << "\nflag 5" << endl;
        Tcritical = T;
        cout << "Tcritical = " << Tcritical << endl;
        flag = 6;
        T = T - 1;
        }

        if(flag==7)
        {
        T = T-0.1; //From flag 6

            for(int l=0; l<n; l++)
            {
            crit_isot << std::fixed << std::setprecision(15) << rho_vec_out[l] << ";" << f_vec_out[l] << ";"
                      << f0_vec_out[l] << ";" << u_vec[l] << ";" << P_vec[l]<< ";" << u_vec_0[l] << ";" << P_vec_0[l] << ";" << T << endl;
            }

        cout << "\nflag 7" << endl;
        T = Tnew+1;

        double b_crit;

            if(exponents==1)
            {
            //double beta_crit = critical_exponents(EdE);
            beta_exponent(&b_crit);
            cout << "beta out of critical = " << b_crit << endl;
            }

        ofstream out_simulation("simulation_record.csv", fstream::app);
        out_simulation << EdE << ";" << cp[0] << ";" << Tcritical << ";" << Pcritical << ";" << rhocritical << ";"
                       << ucritical << ";" << b_crit << ";" << L << ";" << phi_r << endl;
        }

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
