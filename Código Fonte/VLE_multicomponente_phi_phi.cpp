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
#include "mixingrules.h"
#include "EdE.h"
#include "Gibbs.h"
#include "Association.h"

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
VectorXd ln10(nc), logPsat(nc), lnPsat(nc), x(nc), y(nc), phi_liquid_phase(nc), phi_vapor_phase(nc);
VectorXd a0(nc), bCPA(nc), E(nc), beta(nc), c1(nc), X(4*nc), Xl(4*nc), Xv(4*nc);
VectorXd n_v(nc), n_vx(nc), n_vy(nc), yinit(nc);
MatrixXd DELTA (4*nc,4*nc), alfa_NRTL(nc,nc), Aij(nc,nc);
MatrixXd E_row(4*nc,4*nc), E_col(4*nc,4*nc), beta_row(4*nc,4*nc), beta_col(4*nc,4*nc), E_i(4,4), beta_i(4,4), E_iRT(4,4);
VectorXd u_liquid(nc), u_vapor(nc), cond_u(nc);
double tol_u, max_cond_u, u_liquid1, u_vapor1;
int combining_rule, assoc_scheme;
double Tc_a[nc], Pc_a[nc], omega_a[nc], r_a[nc], q_a[nc], q__a[nc], A_a[nc], B_a[nc], C_a[nc]; //Vetores em C++ apenas para gravar valores do arquivo
double E_a[nc], beta_a[nc], a0_a[nc], bCPA_a[nc], c1_a[nc];
double P, T, R, Pnew, sumKx, al, bl, av, bv, Ey;
double tolZv, tolZl, tolSUMKx, tolKx, initialSUMKx, sumKxnew, finalSUMKx, tolX, tolV;
double errorZv, errorZl, errorSUMKx, errorKx, k12, V, dP_dV, rho_l, X1, Vl, Vv, Vt, deltaV;
double Pinit, Vlinit, Vvinit, Zl, Zv, X1l, X1v;
VectorXd ln_Ki(nc), pre_P1(nc), Tr_1(nc), pre_P1_exp(nc), pre_ln_Ki;
double Vl_obj, Vv_obj, Ql, Qv, dP_dVl, dP_dVv;
double log10P, Tb, Tinit, Told;
double G_ex;
VectorXd Tsat(nc), Alog10P(nc), gama(nc), ln_gama(nc);
double init_T, final_T, step, BETCR;
int max_num_iter, counter, stop;
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
double prop[150][18]; //Matrix to hold all data from properties.csv organized

    ifstream file("../Planilhas de análise/properties.csv");

    for(int row = 0; row < 150; ++row)
    {
        string line;
        getline(file, line);
        if ( !file.good() )
            break;
        stringstream iss(line);

        for (int col = 0; col < 18; ++col)
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

//Tolerances
tolZv = 0.0000001; //Erro para convergência de Zv
tolZl = 0.0000001; //Erro para convergência de Zl
tolSUMKx = 0.000001; //Erro para convergência no somatório de Kx
tolKx = 0.00001; //Erro para convergência de Kx
tolX = 0.00000001; //Fraction of non-associating sites tolerance
tolV = 0.0000001; //Volume tolerance

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
    E_i = energy_auto_association(assoc_type[j], E_component, R, T);
    beta_i = volume_auto_association(assoc_type[j], beta_component);

    for (i=0; i<nc; i++)
    {
        E_col.block(4*i,4*j,4,4) << E_i;
        beta_col.block(4*i,4*j,4,4) = beta_i;
    }

}

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
output     << "x1 " << ";" << "y1 " << ";" << "T " << ";" << "P" << ";" << "Vl" << ";"
           << "Vv" << ";" << "sumKx" << ";" << "counter" << ";" << "u_liquid1" << ";" << "u_vapor1" << ";"
           << "X1l" << ";" << "X1v" << ";" << "Zl" << ";" << "Zv" << ";" << "phi_liquid_1"
           << ";" << "phi_liquid_2" << ";" << "phi_vapor_1" << ";" << "phi_vapor_2" << ";"
           << "Vl_obj" << ";" << "Vv_obj" << ";" << "dP/dVl" << ";" << "dP/dVv" << ";" <<
              "G_excess" << endl;

//CÁLCULO PARA COMPONENTE PURO
//============================================================================================================
if(mixture==1)
{

cout << "\nDefine initial Temperature: ";
cin >> init_T;

cout << "\nDefine final Temperature: ";
cin >> final_T;

cout << "\nDefine steps: ";
cin >> step;

T = init_T;
Told = T;

for (T = init_T ; T<=final_T ; T=T+step)
{
 x(0) = 0.999999;
 x(1) = 1-x(0);

E_row = ((E_row.array().log())*Told/T).exp();
E_col = ((E_col.array().log())*Told/T).exp();

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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR);
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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR);
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


Told = T;



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

cout << "T = " << T << endl;
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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vl, Vlinit, a, &Vl_obj, &Ql, &dP_dVl, BETCR);

    phase = 2;
    Xv = volume_function(nc, EdE, phase, y, Xv, EdE_parameters, bv, av, R, T, P, tolV, tolZv, b, combining_rule, beta_row,
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR);
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
                        beta_col, E_row, E_col, alfa, tolX, n_v, &Vv, Vvinit, a, &Vv_obj, &Qv, &dP_dVv, BETCR);
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
