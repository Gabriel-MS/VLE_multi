#ifndef NUMERICAL_H_INCLUDED
#define NUMERICAL_H_INCLUDED
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "interpolation_util.h"
#include <math.h>

#define ITMAX 100
#define EPS 3.0e-8
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//Golden section
#define Rg 0.61803399
#define Cgold (1.0-Rg)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>

//Funtion to return derivative vectors using finite difference method
//Using Eigen type vectors
VectorXd deriv_vec_eigen(VectorXd y, VectorXd x, int deriv)
{
    int SIZE, SIZE2, i;
    SIZE = y.size();
    SIZE2 = SIZE-1;
    VectorXd d(SIZE);

    if(deriv==1)
    {
        d(0) = (y(1)-y(0))/(x(1)-x(0));
        for(i=1; i<SIZE2; i++)
        {
        d(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
        }
        d(SIZE2) = (y(SIZE2)-y(SIZE2-1))/(x(SIZE2)-x(SIZE2-1));
    }

    return d;
}


//Trapezoidal method to calculate area of integral, , uses discrete data points
//Outputs area under function
double trapezoidal_rule(vector<double>& x, vector<double>& y)
{
int i = 0;
int SIZE = x.size();
double sum = 0;
double h = (x[SIZE-1]-x[0])/SIZE;
sum = sum + y[0];
for(i=1; i<SIZE; i++)
{
sum = sum + 2*y[i];
}
sum = sum + y[SIZE-1];
sum = sum * h/2;
return sum;
}


//Trapezoidal method to calculate area of integral, uses discrete data points
//Outputs area under function AND GIVEN INTERVAL
double trapezoidal_interval(vector<double>& x, vector<double>& y, int a, int b)
{
int i = 0;
double sum = 0;
double h = (x[b]-x[a])/(b-a);
sum = sum + y[a];
for(i=a+1; i<b; i++)
{
sum = sum + 2*y[i];
}
sum = sum + y[b];
sum = sum * h/2;
return sum;
}

//Simpson method to calculate area of integral, uses discrete data points
//Outputs area under function
double simpson_rule(vector<double>& x, vector<double>& y)
{
int i = 0;
int SIZE = x.size();
double sum = 0;
double h = (x[SIZE-1]-x[0])/SIZE;
for(i=1; i<SIZE; i++)
{
if(i%2 > 0) sum = sum + 4*y[i];
else sum = sum + 2*y[i];
}
sum = sum + y[0] + y[SIZE-1];
sum = sum * h/3;
return sum;
}

//Simpson method to calculate area of integral, uses discrete data points
//Outputs area under function AND GIVEN INTERVAL
double simpson_interval(vector<double>& x, vector<double>& y, int a, int b)
{
int i = 0;
double sum = 0;
double h = (x[b]-x[a])/(b-a);
for(i=a+1; i<b; i++)
{
if(i%2 > 0) sum = sum + 4*y[i];
else sum = sum + 2*y[i];
}
sum = sum + y[a] + y[b];
sum = sum * h/3;
return sum;
}


//Calculates enclosed area using integral area method, and trapezoidal method
//Enclosed area between two x points and curve, area under curve calculated with trapezoidal method
double enclosed_area_trap(double trap_area, int a, int b)
{
//std::vector <double> x_subv(b-a+1), y_subv(b-a+1); //Cutting subvectors
int i=0;
double Area; //trap_area;
/*
for(i=0; i<b; i++)
{
y_subv[i] = y[a];
x_subv[i] = x[a];
a++;
}
*/
//Area under curve
//trap_area = trapezoidal_rule(x_subv,y_subv);
//Area = trap_area - (y[b]+y[a])/2*(x[b]-x[a]);
return Area;
}


//Calculates enclosed area using integral area method, and trapezoidal method
//Enclosed area between two x points and curve, area under curve calculated with simpson method
double enclosed_area_simp(vector<double>& x, vector<double>& y, int a, int b)
{
std::vector <double> x_subv(b-a), y_subv(b-a); //Cutting subvectors
int i=0;
double Area, trap_area;
for(i=0; i<b; i++)
{
y_subv[i] = y[a];
x_subv[i] = x[a];
a++;
}
//Area under curve
trap_area = simpson_rule(x_subv,y_subv);
Area = trap_area - (y[b]+y[a])/2*(x[b]-x[a]);
return Area;
}


///Finite difference method to calculate approximate first derivative
//x and y vectors must have the same size
//Outputs first derivatives in vector
vector<double> fin_diff_1_vec(vector<double>& x, vector<double>& y)
{
int SIZE, i, SIZE2;
SIZE = x.size();
SIZE2 = SIZE-1;

std::vector<double> y1(SIZE); //Vector to output first derivatives
y1[0] = (y[1]-y[0])/(x[1]-x[0]);
for(i=1; i<SIZE-1; i++)
{
y1[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
}
y1[SIZE2] = (y[SIZE2]-y[SIZE2-1])/(x[SIZE2]-x[SIZE2-1]);

return y1;
}


//Finite difference method to calculate approximate first derivative
//x and y vectors must have the same size
//x MUST have equal spacing between neighbour characters
//Outputs first derivatives in vector
vector<double> fin_diff_2_vec(const vector<double>& x, const vector<double>& y)
{
int SIZE, i;
SIZE = x.size();
std::vector<double> y2(SIZE); //Vector to output first derivatives
y2[0] = (y[2] - 2*y[1] + y[0])/(pow(x[1]-x[0],2)); //Wrong, fix
for(i=1; i<SIZE-1; i++)
{
y2[i] = (y[i+1]-2*y[i]+y[i-1])/(pow(x[i]-x[i-1],2));
}
y2[SIZE-1] = (y[SIZE-1] - 2*y[SIZE-2] + y[SIZE-3])/(pow(x[SIZE-2]-x[SIZE-3],2));
return y2;
}


//Interval method to find root between interval, outputs x value of root
//Using discrete data points, outputs x vector position at root
double bisection_root(vector<double>& x, vector<double>& y, int sup, int inf, double tol)
{
//inf and sup are inferior and superior positions of x vector to contain initial interval
//tol is the tolerance for the method
int cen, counter;
double c;
c = tol+1;
counter = 0;
cout << "sup initial = " << x[sup] << " | inf initial = " << x[inf] << endl;
while(fabs(c)>tol || counter == 100)
{
cen = (sup+inf)/2;
c = y[cen];
if(y[cen]*y[inf] < 0)
{
sup = cen;
}
else
{
inf = cen;
}
if(abs(inf-sup)<2)
{
c = 1e-10;
cout << "inf-sup = " << inf-sup << endl;
}
cout <<"cen = " << cen << " | c = " << c << " | sup = " << x[sup] << " | inf = " << x[inf] << endl;
counter++;
}
return cen;
}


//Cubic spline interpolation, interpolates single values in cubic splines
double cspline(vector<double>& X, vector<double>& Y, double x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = X.size();
double y1; //Output first derivatives
y1 = s(x);
return y1;
}


//Cubic spline interpolation, interpolates vectors in cubic splines
vector<double> cspline_vec(vector<double>& X, vector<double>& Y, vector<double>& x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = x.size();
std::vector<double> y1(SIZE); //Vector to output interpolated values
for(i=1; i<SIZE; i++)
{
y1[i] = s(x[i]);
}
return y1;
}


//Cubic spline derivative, outputs values of 1st derivatives in cubic splines
double cspline_deriv1(vector<double>& X, vector<double>& Y, double x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = X.size();
double y1; //Vector to output first derivatives
y1 = s.deriv(1,x);
return y1;
}


//Cubic spline derivative, outputs vectors of 1st derivatives in cubic splines
vector<double> cspline_deriv1_vec(vector<double>& X, vector<double>& Y, vector<double>& x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = x.size();
std::vector<double> y1(SIZE); //Vector to output first derivatives
for(i=1; i<SIZE; i++)
{
y1[i] = s.deriv(1,x[i]);
}
return y1;
}


//Cubic spline derivative, outputs values of 1st derivatives in cubic splines
double cspline_deriv2(vector<double>& X, vector<double>& Y, double x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = X.size();
double y1; //Vector to output first derivatives
y1 = s.deriv(2,x);
return y1;
}


//Cubic spline derivative, outputs vectors of 1st derivatives in cubic splines
vector<double> cspline_deriv2_vec(vector<double>& X, vector<double>& Y, vector<double>& x)
{
int SIZE, i;
tk::spline s;
s.set_points(X,Y);
SIZE = x.size();
std::vector<double> y1(SIZE); //Vector to output first derivatives
for(i=1; i<SIZE; i++)
{
y1[i] = s.deriv(2,x[i]);
}
return y1;
}

//Calculate phase_coexistence from isotherm and data given by renormalization code
//All vectors must have the same size
double phase_coexistence(double T, vector<double> P, vector<double> u, vector<double> rho, int phase)
{
    int SIZE, i, SIZE2;
    double rho_liq, rho_vap, rho1, drhol, drhov, dPl, dPv, dul, duv, Pliq, Pvap, uliq, uvap, detJ, cond_rhov, cond_rhol, tol;
    SIZE = P.size(); //Measure size of vectors
    SIZE2 = SIZE-1;

    vector<double> dP_drho(SIZE), du_drho(SIZE);
    vector<double> F(2);
    double J[2][2];

    cout << "before cond " << endl;
    tol = 1e-3;
    cond_rhol = tol+1;
    cond_rhov = tol+1;

    //Initial guess for solution
    cout << "before guess" << endl;
    Pliq = P[990];
    Pvap = P[10];
    uliq = u[990];
    uvap = u[10];
    rho_liq = rho[990];
    rho_vap = rho[10];

    //Pressure and chemical potential derivatives for density
    cout << "before deriv" << endl;
    dP_drho[0] = (P[1]-P[0])/(rho[1]-rho[0]);
    for(i=1; i<SIZE-1; i++)
    {
    dP_drho[i] = (P[i+1]-P[i-1])/(rho[i+1]-rho[i-1]);
    }
    dP_drho[SIZE] = (P[SIZE]-P[SIZE2])/(rho[SIZE]-rho[SIZE2]);

    du_drho[0] = (u[1]-u[0])/(rho[1]-rho[0]);
    for(i=1; i<SIZE-1; i++)
    {
    du_drho[i] = (u[i+1]-u[i-1])/(rho[i+1]-rho[i-1]);
    }
    du_drho[SIZE] = (u[SIZE]-u[SIZE2])/(rho[SIZE]-rho[SIZE2]);

    //dP_drho = fin_diff_1_vec(*rho, *P);
    //du_drho = fin_diff_1_vec(*rho, *u);
    cout << "before iter" << endl;
    while(cond_rhol > tol || cond_rhov > tol)
    {
    //Objective function
    F[0] = Pvap - Pliq;
    F[1] = uvap - uliq;

    cout << "before spline" << endl;
    dPl = cspline(rho,dP_drho,rho_liq);
    dPv = cspline(rho,dP_drho,rho_vap);
    dul = cspline(rho,du_drho,rho_liq);
    duv = cspline(rho,du_drho,rho_vap);

    cout << "before J" << endl;
    //jacobian Matrix
    J[0][0] = -dPl;
    J[0][1] = dPv;
    J[1][0] = -dul;
    J[1][1] = duv;
    detJ = J[0][0]*J[1][1] - J[1][0]*J[0][1];

    //Solving J^-1 * F
    drhol = -duv/detJ*F[0]+dPv/detJ*F[1];
    drhov = -dul/detJ*F[0]+dPl/detJ*F[1];

    rho_liq = rho_liq + drhol;
    rho_vap = rho_vap + drhov;

    Pliq = cspline(rho,P,rho_liq);
    Pvap = cspline(rho,P,rho_vap);
    uliq = cspline(rho,u,rho_liq);
    uvap = cspline(rho,u,rho_vap);

    cond_rhol = fabs(drhol/rho_liq);
    cond_rhov = fabs(drhov/rho_vap);
    cout << "cond_rhol = " << cond_rhol << endl;
    }


    if(phase==1)
    {
        rho1 = rho_liq; //outputs liquid density
    }

    else
    {
       rho1 = rho_vap; //outputs vapor density
    }


    return rho1;
}

//Calculate phase_coexistence from isotherm and data given by renormalization code
//All vectors must have the same size
//USING EIGEN TYPE VECTORS
double phase_coexistence_e2(double T, VectorXd P, VectorXd u, VectorXd rho, VectorXd f, int phase)
{
    int SIZE, i, SIZE2;
    double rho_liq, rho_vap, rho1, drhol, drhov, dPl, dPv, dul, duv, Pliq, Pvap, uliq, uvap, detJ, cond_rhov, cond_rhol, tol;
    double fliq, fvap;
    SIZE = P.size(); //Measure size of vectors
    SIZE2 = SIZE-1;

    //std::vector<double> *rho2, *dP_drho2, *du_drho2, *P2, *u2;
    VectorXd dP_drho(SIZE), du_drho(SIZE);
    VectorXd F(2);
    MatrixXd J(2,2);

    //cout << "before cond " << endl;
    tol = 1e-4;
    cond_rhol = tol+1;
    cond_rhov = tol+1;

    //Initial guess for solution
    //cout << "before guess" << endl;
    Pliq = P(700);
    Pvap = P(50);
    uliq = u(700);
    uvap = u(50);
    rho_liq = rho(700);
    rho_vap = rho(50);

    //Pressure and chemical potential derivatives for density
    //cout << "before deriv" << endl;

    int one;
    //cout << "size = " << SIZE << endl;
    dP_drho = deriv_vec_eigen(P, rho, 1);
    du_drho = deriv_vec_eigen(u, rho, 1);
    //cout << "before iter " << endl;
    //cin >> one;

/*
    for(i=0;i<1000;i++)
    {
        cout << "i = " << i << " / " << dP_drho(i) - (P(i+1)-P(i-1))/(rho(i+1)-rho(i-1)) << endl;
    }
*/

    std::vector<double> dP_drho2(1000), du_drho2(1000), rho2(1000), P2(1000), u2(1000), f2(1000);
    //cout << "initial guess = " << rho_liq << " / " << rho_vap << endl;

for(i=0;i<1000;i++)
{
    dP_drho2[i] = dP_drho(i);
    du_drho2[i] = du_drho(i);
    rho2[i] = rho(i);
    P2[i] = P(i);
    u2[i] = u(i);
    f2[i] = f(i);
}

//cout << "before conditions" << endl;
//cout << "liq / vap = " << rho_liq << " / " << rho_vap << endl;
    i = 0;
    while((cond_rhol > tol || cond_rhov > tol) && (i<500))
    {
    //Objective function
    F(0) = Pvap - Pliq;
    F(1) = uvap - uliq;

    //cout << "before spline" << endl;
    dPl = cspline(rho2,dP_drho2,rho_liq);
    dPv = cspline(rho2,dP_drho2,rho_vap);
    dul = cspline(rho2,du_drho2,rho_liq);
    duv = cspline(rho2,du_drho2,rho_vap);

    //cout << "before J" << endl;
    //jacobian Matrix
    J(0,0) = -dPl;
    J(0,1) = dPv;
    J(1,0) = -dul;
    J(1,1) = duv;
    detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    //cout << "J = " << J << endl;
    //cout << "detJ = " << detJ << endl;

    //Solving J^-1 * F
    drhol = -duv/detJ*F(0)+dPv/detJ*F(1);
    drhov = -dul/detJ*F(0)+dPl/detJ*F(1);

    if(i<2)
    {
    rho_liq = rho_liq + drhol/50;
    rho_vap = rho_vap + drhov/50;
    }

    else
    {
    rho_liq = rho_liq + drhol/10;
    rho_vap = rho_vap + drhov/10;
    }

    fliq = cspline(rho2,f2,rho_liq);
    fvap = cspline(rho2,f2,rho_vap);
    uliq = cspline_deriv1(rho2,f2,rho_liq);
    uvap = cspline_deriv1(rho2,f2,rho_vap);
    Pliq = -fliq + rho_liq*uliq;
    Pvap = -fvap + rho_vap*uvap;

    //Pliq = cspline(rho2,P2,rho_liq);
    //Pvap = cspline(rho2,P2,rho_vap);
    //uliq = cspline(rho2,u2,rho_liq);
    //uvap = cspline(rho2,u2,rho_vap);

    //cond_rhol = fabs(drhol/rho_liq);
    //cond_rhov = fabs(drhov/rho_vap);
    //cond_rhol = fabs(drhol);
    //cond_rhov = fabs(drhov);
    cond_rhol = fabs(F(0));
    cond_rhov = fabs(F(1));
    //cout << "cond_rhol = " << cond_rhol << " " << cond_rhov << endl;
    cout << "liq / vap = " << rho_liq << " / " << rho_vap << " " << i << endl;
    //cin >> one;
    i = i + 1;
    }


    if(phase==1)
    {
        rho1 = rho_liq; //outputs liquid density
    }

    else
    {
       rho1 = rho_vap; //outputs vapor density
    }

cout << "rho1 = " << rho1 << endl;
    return rho1;
}

//Calculate phase_coexistence from isotherm and data given by renormalization code
//All vectors must have the same size
//USING EIGEN TYPE VECTORS
double phase_coexistence_e3(double T, VectorXd P, VectorXd u, VectorXd rho, VectorXd f, int posl, int posv, int phase)
{
    int SIZE, i, SIZE2;
    double rho_liq, rho_vap, rho1, drhol, drhov, dPl, dPv, dul, duv, Pliq, Pvap, uliq, uvap, detJ, cond_rhov, cond_rhol, tol;
    double fliq, fvap;
    SIZE = P.size(); //Measure size of vectors
    SIZE2 = SIZE-1;

    //std::vector<double> *rho2, *dP_drho2, *du_drho2, *P2, *u2;
    VectorXd dP_drho(SIZE), du_drho(SIZE);
    VectorXd F(2);
    MatrixXd J(2,2);

    //cout << "before cond " << endl;
    tol = 1e-4;
    cond_rhol = tol+1;
    cond_rhov = tol+1;

    //Initial guess for solution
    //cout << "before guess" << endl;
    Pliq = P(posl);
    Pvap = P(posv);
    uliq = u(posl);
    uvap = u(posv);
    rho_liq = rho(posl);
    rho_vap = rho(posv);

    //Pressure and chemical potential derivatives for density
    //cout << "before deriv" << endl;

    int one;
    //cout << "size = " << SIZE << endl;
    dP_drho = deriv_vec_eigen(P, rho, 1);
    du_drho = deriv_vec_eigen(u, rho, 1);
    //cout << "before iter " << endl;
    //cin >> one;

/*
    for(i=0;i<1000;i++)
    {
        cout << "i = " << i << " / " << dP_drho(i) - (P(i+1)-P(i-1))/(rho(i+1)-rho(i-1)) << endl;
    }
*/

    std::vector<double> dP_drho2(1000), du_drho2(1000), rho2(1000), P2(1000), u2(1000), f2(1000);
    //cout << "initial guess = " << rho_liq << " / " << rho_vap << endl;

for(i=0;i<1000;i++)
{
    dP_drho2[i] = dP_drho(i);
    du_drho2[i] = du_drho(i);
    rho2[i] = rho(i);
    P2[i] = P(i);
    u2[i] = u(i);
    f2[i] = f(i);
}

//cout << "before conditions" << endl;
//cout << "liq / vap = " << rho_liq << " / " << rho_vap << endl;
    i = 0;
    while((cond_rhol > tol || cond_rhov > tol) && (i<500))
    {
    //Objective function
    F(0) = Pvap - Pliq;
    F(1) = uvap - uliq;

    //cout << "before spline" << endl;
    dPl = cspline(rho2,dP_drho2,rho_liq);
    dPv = cspline(rho2,dP_drho2,rho_vap);
    dul = cspline(rho2,du_drho2,rho_liq);
    duv = cspline(rho2,du_drho2,rho_vap);

    //cout << "before J" << endl;
    //jacobian Matrix
    J(0,0) = -dPl;
    J(0,1) = dPv;
    J(1,0) = -dul;
    J(1,1) = duv;
    detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    //cout << "J = " << J << endl;
    //cout << "detJ = " << detJ << endl;

    //Solving J^-1 * F
    drhol = -duv/detJ*F(0)+dPv/detJ*F(1);
    drhov = -dul/detJ*F(0)+dPl/detJ*F(1);

    if(i<2)
    {
    rho_liq = rho_liq + drhol/50;
    rho_vap = rho_vap + drhov/50;
    }

    else
    {
    rho_liq = rho_liq + drhol/5; //ERA 10
    rho_vap = rho_vap + drhov/5; //ERA 10
    }

    fliq = cspline(rho2,f2,rho_liq);
    fvap = cspline(rho2,f2,rho_vap);
    uliq = cspline_deriv1(rho2,f2,rho_liq);
    uvap = cspline_deriv1(rho2,f2,rho_vap);
    Pliq = -fliq + rho_liq*uliq;
    Pvap = -fvap + rho_vap*uvap;

    //Pliq = cspline(rho2,P2,rho_liq);
    //Pvap = cspline(rho2,P2,rho_vap);
    //uliq = cspline(rho2,u2,rho_liq);
    //uvap = cspline(rho2,u2,rho_vap);

    //cond_rhol = fabs(drhol/rho_liq);
    //cond_rhov = fabs(drhov/rho_vap);
    //cond_rhol = fabs(drhol);
    //cond_rhov = fabs(drhov);
    cond_rhol = fabs(F(0));
    cond_rhov = fabs(F(1));
    //cout << "cond_rhol = " << cond_rhol << " " << cond_rhov << endl;
    cout << "liq / vap = " << rho_liq << " / " << rho_vap << " " << i << endl;
    //cin >> one;
    i = i + 1;
    }


    if(phase==1)
    {
        rho1 = rho_liq; //outputs liquid density
    }

    else
    {
       rho1 = rho_vap; //outputs vapor density
    }

cout << "rho1 = " << rho1 << endl;
    return rho1;
}

//Calculate phase_coexistence from isotherm and data given by renormalization code
//All vectors must have the same size
//USING EIGEN TYPE VECTORS - integrate area method with newton-raphson
double phase_coexistence_e(double T, VectorXd P, VectorXd f, VectorXd rho, int phase)
{
    int SIZE, i, SIZE2;
    SIZE = P.size(); //Measure size of vectors
    SIZE = 1000;
    SIZE2 = SIZE-1;

    double V1, V1p, V1m, V2, V2p, V2m, dAdV1, dAdV2, dAdV1m, dAdV1p, dAdV2m, dAdV2p, Vout, tol;
    double condV1, condV2, cond_rhol, cond_rhov, AV1, AV1m, AV1p, AV2, AV2m, AV2p;
    double detJ, P1, P2, dV1, dV2;
    std::vector<double> Vn(1000), An(1000), dAdVn(1000), Pn(1000);
    VectorXd A(1000), dAdV(1000), V(1000);
    VectorXd F(2), dV(2);
    MatrixXd J(2,2);

    dAdV = deriv_vec_eigen(A, V, 1);

    for(i=0;i<1000;i++)
    {
        A(i) = f(i)/rho(i);
        V(i) = 1/rho(i);
        Vn[i] = V(i);
        An[i] = A(i);
        dAdVn[i] = dAdV(i);
        Pn[i] = P(i);
    }


    cout << "before cond " << endl;
    tol = 1e-3;
    condV1 = tol+1;
    condV2 = tol+1;

    //Initial guess for solution
    cout << "before guess" << endl;
    V2 = V(10);
    V1 = V(990);
    dAdV2 = dAdV(10);
    dAdV1 = dAdV(990);

    condV1 = tol+1;
    condV2 = tol+1;

    //Pressure and chemical potential derivatives for density
    cout << "before deriv" << endl;

    int one;
    cout << "size = " << SIZE << endl;

    cout << "before iter " << endl;


cout << "before conditions" << endl;

    while(condV1 > tol || condV2 > tol)
    {
     dAdV1 = cspline(Vn, dAdVn, V1);
     dAdV1m = cspline(Vn, dAdVn, V1-1e-4);
     dAdV1p = cspline(Vn, dAdVn, V1+1e-4);
     dAdV2 = cspline(Vn, dAdVn, V2);
     dAdV2m = cspline(Vn, dAdVn, V2-1e-4);
     dAdV2p = cspline(Vn, dAdVn, V2+1e-4);

     AV1 = cspline(Vn, An, V1);
     AV1m = cspline(Vn, An, V1-1e-4);
     AV1p = cspline(Vn, An, V1+1e-4);
     AV2 = cspline(Vn, An, V1);
     AV2m = cspline(Vn, An, V1-1e-4);
     AV2p = cspline(Vn, An, V1+1e-4);

    //Objective function
    F(0) = (V2-V1)/2*dAdV1 + (AV2-AV1)/2;
    F(1) = (V2-V1)/2*dAdV2 + (AV2-AV1)/2;

    cout << "before J" << endl;
    //jacobian Matrix
    J(0,0) = ((V2-V1+1e-4)/2*dAdV1p + (AV2-AV1p)/2  -  (V2-V1-1e-4)/2*dAdV1m + (AV2-AV1m)/2)/(2e-4);
    J(0,1) = ((V2-V1+1e-4)/2*dAdV1 + (AV2p-AV1)/2  -  (V2-V1-1e-4)/2*dAdV1 + (AV2m-AV1)/2)/(2e-4);
    J(1,0) = ((V2-V1+1e-4)/2*dAdV2 + (AV2-AV1p)/2  -  (V2-V1-1e-4)/2*dAdV2 + (AV2-AV1m)/2)/(2e-4);
    J(1,1) = ((V2-V1+1e-4)/2*dAdV2p + (AV2p-AV1p)/2  -  (V2-V1-1e-4)/2*dAdV2m + (AV2-AV1m)/2)/(2e-4);
    detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    cout << "J = " << J << endl;
    cout << "detJ = " << detJ << endl;

    //Solving J^-1 * F
    dV = J.inverse()*F;

    V1 = V1 + dV(0);
    V2 = V2 + dV(1);

    P1 = cspline(Vn, Pn, V1);
    P2 = cspline(Vn, Pn, V2);

    condV1 = fabs(dV1/V1);
    condV2 = fabs(dV2/V2);
    cout << "cond_rhol = " << condV1 << " " << condV2 << endl;
    cout << "liq / vap = " << 1/V1 << " / " << 1/V2 << endl;
    cin >> one;
    }


    if(phase==1)
    {
        Vout = V1; //outputs liquid density
    }

    else
    {
       Vout = V2; //outputs vapor density
    }

cout << "V out = " << 1/Vout << endl;
    return Vout;
}

//Given function and interval, find root using brent method
//Numerical Recipes in C (R)
double zbrent(double (*func)(double), double x1, double x2, double tol)
{
int iter;
double a=x1,b=x2,c=x2,d,e,min1,min2;
double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
fc=fb;
for (iter=1;iter<=ITMAX;iter++) {
if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
c=a;
fc=fa;
e=d=b-a;
}
if (fabs(fc) < fabs(fb)) {
a=b;
b=c;
c=a;
fa=fb;
fb=fc;
fc=fa;
}
tol1=2.0*EPS*fabs(b)+0.5*tol;
xm=0.5*(c-b);
if (fabs(xm) <= tol1 || fb == 0.0) return b;
if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
s=fb/fa;
if (a == c) {
p=2.0*xm*s;
q=1.0-s;
} else {
q=fa/fc;
r=fb/fc;
p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
q=(q-1.0)*(r-1.0)*(s-1.0);
}
if (p > 0.0) q = -q;
p=fabs(p);
min1=3.0*xm*q-fabs(tol1*q);
min2=fabs(e*q);
if (2.0*p < (min1 < min2 ? min1 : min2)) {
e=d;
d=p/q;
} else {
d=xm;
e=d;
}
} else {
d=xm;
e=d;
}
a=b;
fa=fb;
if (fabs(d) > tol1)
b+=d;
else
b += SIGN(tol1,xm);
fb=(*func)(b);
}
return 0.0;
}

//Given function and interval, find root using brent method
//Modified to fit Maxwell construction function
//Numerical Recipes in C (R)
double zbrentm(double (*func)(vector<double>&, vector<double>&, double), double x1, double x2,
               double tol, vector<double>& v1, vector<double>& v2)
{
int iter;
double a=x1,b=x2,c=x2,d,e,min1,min2;
double fa=(*func)(v1,v2,a),fb=(*func)(v1,v2,b),fc,p,q,r,s,tol1,xm;
if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
fc=fb;
for (iter=1;iter<=ITMAX;iter++) {
if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
c=a;
fc=fa;
e=d=b-a;
}
if (fabs(fc) < fabs(fb)) {
a=b;
b=c;
c=a;
fa=fb;
fb=fc;
fc=fa;
}
tol1=2.0*EPS*fabs(b)+0.5*tol;
xm=0.5*(c-b);
if (fabs(xm) <= tol1 || fb == 0.0) return b;
if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
s=fb/fa;
if (a == c) {
p=2.0*xm*s;
q=1.0-s;
} else {
q=fa/fc;
r=fb/fc;
p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
q=(q-1.0)*(r-1.0)*(s-1.0);
}
if (p > 0.0) q = -q;
p=fabs(p);
min1=3.0*xm*q-fabs(tol1*q);
min2=fabs(e*q);
if (2.0*p < (min1 < min2 ? min1 : min2)) {
e=d;
d=p/q;
} else {
d=xm;
e=d;
}
} else {
d=xm;
e=d;
}
a=b;
fa=fb;
if (fabs(d) > tol1)
b+=d;
else
b += SIGN(tol1,xm);
fb=(*func)(v1,v2,b);
}
return 0.0;
}

//Finds maximum integer position of dPdV inside binodal curve
int bin_max(vector<double>& dP_dV)
{
double pmax_cond;
int i=50;
do
{
 pmax_cond = dP_dV[i];
 i++;
//cout << i << endl;
}while(pmax_cond>0);
    return i-1;
}

//Finds minimum integer position of dPdV inside binodal curve
int bin_min(vector<double>& dP_dV)
{
int j=900;
double pmin_cond;
do
{
 pmin_cond = dP_dV[j];
 j--;
//cout << j << endl;
//cout << "j = " << j << " / pmin = " << pmin_cond << " / rho = " << rho_vec_out[j] << endl;
}while(pmin_cond>0);
    return j+1;
}

//Finds maximum integer position of dPdV inside binodal curve
int bin_min_seed(vector<double>& dP_dV, int seed)
{
double pmin_cond;
int j=seed+1;
do
{
 pmin_cond = dP_dV[j];
 j++;
}while(pmin_cond<0);
    return j-1;

}

//Uses Regula Falsi method with cubic spline interpolation to find root in given interval
double falsi_spline(vector<double>& x, vector<double>& y, double a, double b, double tol)
{
    double ya, yb, yc, c;
    int i;
    yc = tol+1;

    while( fabs(yc)>tol )
    {
    ya = cspline(x, y, a);
    yb = cspline(x, y, b);
    c = b - yb*(a-b)/(ya-yb);
    yc = cspline(x, y, c);
    //cout << "falsi ya yb yc = " << ya << " / " << yb << " / " << yc << endl;
    //cout << "falsi a b c = " << a << " / " << b << " / " << c << endl;
    if(ya*yc<0) b = c;
    else a = c;
    //cout << "AFTER falsi a b c = " << a << " / " << b << " / " << c << endl;
    //cin >> i;
    }

    return c;
}

//Uses Secant method with cubic spline interpolation to find root in given interval
double secant_spline(vector<double>& x, vector<double>& y, double a, double b, double tol)
{
    double c, yc, c_new, yc_new, c_old, yc_old;
    int i;
    yc = tol+1;

    c_old = a;
    yc_old = cspline(x,y,a);
    c = (a+b)/2;
    yc = cspline(x,y,c);

    while( fabs(yc_new)>tol )
    {
    c_new = c - yc*(c_old-c)/(yc_old-yc);
    yc_new = cspline(x,y,yc_new);

    c = c_new;
    yc = yc_new;
    c_old = c;
    yc_old = yc;
    }

    return c;
}

//Uses Secant method with cubic spline interpolation to find root in given interval
double newton_spline(vector<double>& x, vector<double>& y, double c, double ans, double tol)
{
    double F, dF, Fold, yc;
    int i;
    F = tol+1;
    dF = 1000;

    while( fabs(F)>tol )
    {
        Fold = F;
        c = c - F/dF;
        yc = cspline(x,y,c);
        F = fabs(yc-ans);
        dF = F-Fold;
    }

    return c;
}


//Calculates Area of given isotherm (helmholtz energy density f, density rho and chemical potential u)
//All inputs should be provided dimensionless
//Uses Cubic Spline to interpolate data points
double maxwell_cons(vector<double>& rho, vector<double>& P, double P1)
{
    std::vector<double> dPdV(1000), rho_a1(1000), rho_a2(1000), P_a1(1000), P_a2(1000), Pf(1000);
    int SIZE, i;
    double Area1, Area2, Area;
    double max1, min1, rho_max, rho_min, P_max, P_min, rho1, rho2, rho3, rhoz1, rhoz2, h1, h2;

    SIZE = rho.size();

    //Calculating vector with first derivatives
    dPdV = fin_diff_1_vec(rho, P);

    //Find max and min pressure of isotherm inside binodal curve
    max1 = bin_max(dPdV);
    min1 = bin_min(dPdV);
    rho_max = rho[max1];
    rho_min = rho[min1];
    P_max = P[max1];
    P_min = P[min1];

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(P_min<0)
        {
            Pf[i] = P[i] + P_min - P1;
        }

        else
        {
            Pf[i] = P[i] - P1;
        }
    }

    //Find density roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rho_max, 1e-3);
    rho2 = falsi_spline(rho, Pf, rho_max, rho_min, 1e-3);
    rho3 = falsi_spline(rho, Pf, rho_min, rho[700], 1e-3);

    //Find new values of pressure and density between rho1 and rho3
    rhoz1 = rho1;
    rhoz2 = rho2;
    h1 = (rho2-rho1)/1000;
    h2 = (rho3-rho2)/1000;

    rho_a1[0] = rho1;
    rho_a1[999] = rho2;
    rho_a2[0] = rho2;
    rho_a2[999] = rho3;

    P_a1[0] = P1;
    P_a1[999] = P1;
    P_a2[0] = P1;
    P_a2[999] = P1;

    for(i=1; i<999; i++)
    {
    rho_a1[i] = rhoz1;
    rho_a2[i] = rhoz2;
    rhoz1 = rhoz1 + h1;
    rhoz2 = rhoz2 + h2;
    }

    for(i=1; i<999; i++)
    {
    P_a1[i] = cspline(rho, P, rho_a1[i]);
    P_a2[i] = cspline(rho, P, rho_a2[i]);
    }

    //Calculates Maxwell area
    Area1 = simpson_rule(rho_a1, P_a1);
    Area1 = Area1 - P1*(rho2-rho1);

    Area2 = P1*(rho3-rho2);
    Area2 = Area2 - simpson_rule(rho_a2, P_a2);

    Area = Area1 - Area2;
    Area = fabs(Area);

    return Area;
}

//Calculates Area of given isotherm (helmholtz energy density f, density rho and chemical potential u)
//All inputs should be provided dimensionless
//Uses Cubic Spline to interpolate data points
double maxwell_cons2(vector<double>& rho, vector<double>& P, double P1)
{
    std::vector<double> dPdV(1000), rho_a1(1000), rho_a2(1000), P_a1(1000), P_a2(1000), Pf(1000);
    int SIZE, i;
    double Area1, Area2, Area;
    double max1, min1, rho_max, rho_min, P_max, P_min, rho1, rho2, rho3, rhoz1, rhoz2, h1, h2;

    SIZE = rho.size();

    //Calculating vector with first derivatives
    dPdV = fin_diff_1_vec(rho, P);

    //Find max and min pressure of isotherm inside binodal curve
    max1 = bin_max(dPdV);
    min1 = bin_min(dPdV);
    rho_max = rho[max1];
    rho_min = rho[min1];
    P_max = P[max1];
    P_min = P[min1];

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(P_min<0)
        {
            Pf[i] = P[i] + P_min - P1;
        }

        else
        {
            Pf[i] = P[i] - P1;
        }
    }

    //Find density roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rho_max, 1e-3);
    rho2 = falsi_spline(rho, Pf, rho_max, rho_min, 1e-3);
    rho3 = falsi_spline(rho, Pf, rho_min, rho[700], 1e-3);

    //Find new values of pressure and density between rho1 and rho3
    rhoz1 = rho1;
    rhoz2 = rho2;
    h1 = (rho3-rho1)/1000;

    rho_a1[0] = rho1;
    rho_a1[999] = rho3;

    P_a1[0] = P1;
    P_a1[999] = P1;

    for(i=1; i<999; i++)
    {
    rho_a1[i] = rhoz1;
    rhoz1 = rhoz1 + h1;
    }

    for(i=1; i<999; i++)
    {
    P_a1[i] = cspline(rho, P, rho_a1[i]);
    }

    //Calculates Maxwell area
    Area1 = simpson_rule(rho_a1, P_a1);
    Area = Area1 - P1*(rho3-rho1);
    Area = fabs(Area);

    //cout << "Area = " << Area << endl;

    return Area;
}


double area_helm(vector<double>& V, vector<double>& A, vector<double>& P, double P0)
{
    double Area, A1, A2, V1, V2, Vmax, Vmin, Pmax, Pmin, Vh, Vz, Az, h;
    double rhomax, rhomin, rho1, rho2, rhoh, rhoz;
    std::vector<double> Pf(1000), Vp(1000), Ap(1000), dPdV(1000), rho(1000), rhop(1000), dPdrho(1000);
    int i, SIZE, max1, min1;

    SIZE = V.size();
    cout << "area in" << endl;

    for(i=0; i<SIZE; i++)
    {
        rho[i] = 1/V[i];
    }

    //Calculating vector with first derivatives
    dPdrho = fin_diff_1_vec(rho, P);

    //Find max and min pressure of isotherm inside binodal curve
    max1 = bin_max(dPdrho);
    min1 = bin_min(dPdrho);
    rhomax = rho[max1];
    rhomin = rho[min1];
    Pmax = P[max1];
    Pmin = P[min1];

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(Pmin<0)
        {
            Pf[i] = P[i] + Pmin - P0;
        }

        else
        {
            Pf[i] = P[i] - P0;
        }
    }

    cout << "before falsi" << endl;
    //Find volume roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rhomax, 1e-3);
    V1 = 1/rho1;
    rho2 = falsi_spline(rho, Pf, rhomin, rho[700], 1e-3);
    V2 = 1/rho2;
    A1 = cspline(rho, A, rho1);
    A2 = cspline(rho, A, rho2);

    //Find new values of pressure and density between rho1 and rho3
    rhoh = rho1;
    h = (rho2-rho1)/100;

    rhop[0] = rho1;
    rhop[99] = rho2;

    Vp[0] = V1;
    Vp[99] = V2;

    Ap[0] = A1;
    Ap[99] = A2;

    for(i=1; i<99; i++)
    {
    rhop[i] = rhoz;
    Vp[i] = 1/rhop[i];
    Ap[i] = cspline(rho, A, rhop[i]);
    rhoz = rhoz + h;
    //cout << "p = " << Vp[i] << Ap[i] << endl;
    }

    //Calculate Area between V1 and V2
    Area = simpson_rule(Vp,Ap);

    cout << "Area = " << Area << endl;
    cout << "A1 = " << A1 << endl;
    cout << "A2 = " << A2 << endl;
    cout << "V1 = " << V1 << endl;
    cout << "V2 = " << V2 << endl;
    cout << "rho1 = " << rho1 << endl;
    cout << "rho2 = " << rho2 << endl;
    cout << "rhomax = " << rhomax << endl;
    cout << "rhomin = " << rhomin << endl;

    Area = Area - (A2 + A1)/2*(V1-V2); //V1 and V2 switched because V1 is bigger

    cout << "area out = " << Area << endl;
    cin >> i;

    return Area;
}

//golden section search to find maximum value
double golden(double ax, double bx, double cx, double (*f)(vector<double>&, vector<double>&, vector<double>&,
              double), double tol, double *xmin, vector<double>& v1, vector<double>& v2, vector<double>& v3)
{
    double f1,f2,x0,x1,x2,x3;

    x0=ax;
    x3=cx;

    if (fabs(cx-bx) > fabs(bx-ax))
    {
    x1=bx;
    x2=bx+Cgold*(cx-bx);
    }

    else
    {
    x2=bx;
    x1=bx-Cgold*(bx-ax);
    }

    f1=(*f)(v1,v2,v3,x1);
    f2=(*f)(v1,v2,v3,x2);

    while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)))
    {
        if (f2 > f1) //Original numerical recipes code was " < ", for minimum search
        {
        SHFT3(x0,x1,x2,Rg*x1+Cgold*x3)
        SHFT2(f1,f2,(*f)(v1,v2,v3,x2))
        }

        else
        {
        SHFT3(x3,x2,x1,Rg*x2+Cgold*x0)
        SHFT2(f2,f1,(*f)(v1,v2,v3,x1))
        }
    }

    if (f1 > f2)  //Original numerical recipes code was " < ", for minimum search
    {
    *xmin=x1;
    return f1;
    }

    else
    {
    *xmin=x2;
    return f2;
    }
}

//Calcualtes phase coexistence densities and pressure
//Using area method, searching maximum area with golden section
vector<double> dens_area(vector<double>& V, vector<double>& A, vector<double>& P)
{
    int i, SIZE;
    std::vector<double> dens(3), dPdV(1000), Pf(1000), rho(1000);
    double V1, V2, Vmax, Vmin, Pmax, Pmin, P0, P1, P2, Pg, *xmin;
    int max1, min1;

    SIZE = V.size();

    for(i=0; i<SIZE; i++)
    {
        rho[i] = 1/V[i];
    }

    dPdV = fin_diff_1_vec(rho, P);

    //Find max and min of isotherm inside binodal region
    max1 = bin_max(dPdV);
    min1 = bin_min(dPdV);
    Vmax = V[max1];
    Vmin = V[min1];
    Pmax = P[max1];
    Pmin = P[min1];

    //Initial guess for pressure
    P0 = (Pmax + Pmin)/2;
    P1 = Pmax; //range of pressure
    P2 = Pmin; //range of pressure

   //Brent method to use Maxwell Construction Method
    Pg = golden(P1*0.95, P0, P2*1.05, area_helm, 1e-5, xmin, V, A, P);
    Pg = *xmin;
    //cout << "Pz = " << Pz << endl;
    cout << "after golden" << endl;

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(Pmin<0)
        {
            Pf[i] = P[i] + Pmin - Pg;
        }

        else
        {
            Pf[i] = P[i] - Pg;
        }
    }

    //Find initial volume guesses
    V1 = falsi_spline(V, Pf, V[0], Vmax, 1e-3);
    V2 = falsi_spline(V, Pf, Vmin, V[700], 1e-3);

    dens[0] = 1/V1;
    dens[1] = 1/V2;
    dens[3] = Pg;

    return dens;
}

//Function to c3culate density of coexistent phases using maxwell construction
//Input is dimensionless isotherm
vector<double> dens_maxwell(vector<double>& rho, vector<double>& P, vector<double>& f, double tol)
{
    int i, SIZE;
    std::vector<double> dens(7);
    std::vector<double> dPdV(1000), Pf(1000);
    double f1, f2;
    double rho1, rho2, rho_max, rho_min, P_max, P_min, P0, P1, P2, Pz, u1, u2;
    int max1, min1;

    SIZE = rho.size();
    dPdV = fin_diff_1_vec(rho, P);

    //Find max and min of isotherm inside binodal region
    max1 = bin_max(dPdV);
    min1 = bin_min(dPdV);
    rho_max = rho[max1];
    rho_min = rho[min1];
    P_max = P[max1];
    P_min = P[min1];

    //Initial guess for pressure and density
    P0 = (P_max + P_min)/2;
    P1 = P_max; //range of pressure
    P2 = P_min; //range of pressure

    if (P_min>P_max)
    {
        min1 = bin_min_seed(dPdV,max1);
        rho_min = rho[min1];
        P_min = P[min1];
    }

    if(P2<0) P2 = 1e-4;

    //cout << "before brent" << P1*0.99999 << " / " << P2*1.00001 << endl;

    //Brent method to use Maxwell Construction Method
    if(P2>1e-3)
    {
    Pz = zbrentm(maxwell_cons2, P1*0.99999, P2*1.00001, tol, rho, P);

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(P_min<0)
        {
            Pf[i] = P[i] + P_min - Pz;
        }

        else
        {
            Pf[i] = P[i] - Pz;
        }
    }


    //Find density roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rho_max, 1e-3);
    rho2 = falsi_spline(rho, Pf, rho_min, rho[700], 1e-3);

    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;
    }

    //output vector of densities
    dens[0] = rho1;
    dens[1] = rho2;
    dens[2] = Pz;
    dens[3] = P1;
    dens[4] = P2;
    dens[5] = u1;
    dens[6] = u2;
    return dens;
}

//Function to calculate density of coexistent phases using maxwell construction
//Input is dimensionless isotherm
vector<double> dens_maxwell_seed(vector<double>& rho, vector<double>& P, double tol, vector<double>& seed)
{
    int i, SIZE;
    std::vector<double> dens(3);
    std::vector<double> dPdV(1000), Pf(1000);
    double rho1, rho2, rho_max, rho_min, P_max, P_min, P0, P1, P2, Pz, u1, u2;
    int max1, min1;

    SIZE = rho.size();
    dPdV = fin_diff_1_vec(rho, P);

    //Find max and min of isotherm inside binodal region
    max1 = bin_max(dPdV);
    min1 = bin_min(dPdV);
    rho_max = rho[max1];
    rho_min = rho[min1];
    P_max = P[max1];
    P_min = P[min1];

    //Initial guess for pressure and density
    P0 = (P_max + P_min)/2;
    P1 = P_max; //range of pressure
    P2 = P_min; //range of pressure

    //Brent method to use Maxwell Construction Method
    Pz = zbrentm(maxwell_cons2, P1*0.99999, P2*1.00001, tol, rho, P);

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE; i++)
    {
        if(P_min<0)
        {
            Pf[i] = P[i] + P_min - Pz;
        }

        else
        {
            Pf[i] = P[i] - Pz;
        }
    }


    //Find density roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rho_max, 1e-3);
    rho2 = falsi_spline(rho, Pf, rho_min, rho[600], 1e-3);

    //output vector of densities
    dens[0] = rho1;
    dens[1] = rho2;
    dens[2] = Pz;
    return dens;
}

//Function to calculate phase coexistence densities
//Using Newton method
vector<double> dens_newt(vector<double>& rho, vector<double>& f, vector<double>& P, vector<double>& u, double tol)
{
    std::vector<double> dPdrho(1000), Pf1(1000), Pf2(1000), rhov(6);
    int i, max1, min1;
    double rhomax, rhomin, Pmax, Pmin, P1, P2, f1, f2, u1, u2, du1, du2, dP1, dP2, detJ, drho1, drho2;
    double rho1, rho2, area;

    //Calculating vector with first derivatives
    dPdrho = fin_diff_1_vec(rho,P);

    //Find max and min pressure of isotherm inside binodal curve
    max1 = bin_max(dPdrho);
    min1 = bin_min(dPdrho);
    rhomax = rho[max1];
    rhomin = rho[min1];
    Pmax = P[max1];
    Pmin = P[min1];

    if (Pmin>Pmax)
    {
        cout << "entered min seed " << Pmin << " / " << Pmax << endl;
        min1 = bin_min_seed(dPdrho,max1);
        rhomin = rho[min1];
        Pmin = P[min1];
    }

    if (Pmin<0) Pmin=1e-3;

    for(i=0; i<1000; i++)
    {
        if(Pmin<0) Pf1[i] = P[i] + Pmin;
        else Pf1[i] = P[i] - Pmin;

        Pf2[i] = P[i] - Pmax;
    }

    //cout << "Pmax/min = " << Pmax << " / " << Pmin << endl;
    //cout << "max/min = " << max1 << " / " << min1 << endl;
    //cout << "rho = " << rhomax << " / " << rhomin << endl;
    //Find initial guess for densities
    rho1 = falsi_spline(rho, Pf1, rho[0], rhomax, 1e-3);
    //cout << "rho1 = " << rho1 << endl;
    rho2 = falsi_spline(rho, Pf2, rhomin, rho[700], 1e-3);
    //cout << "rho2 = " << rho2 << endl;

    //Solve newton-raphson system
    drho1 = tol+1;
    int counter = 0;
    while((fabs(drho1)>tol || fabs(drho2)>tol) && (fabs(Pmax-Pmin)>1e-3))
    {
    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;

    du1 = cspline_deriv1(rho,u,rho1);
    du2 = cspline_deriv1(rho,u,rho2);
    dP1 = cspline_deriv1(rho,P,rho1);
    dP2 = cspline_deriv1(rho,P,rho2);
    detJ = -dP2*du1+dP1*du2;

    drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2);
    drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2);

    rho1 = rho1 + drho1/10;
    rho2 = rho2 + drho2/10;
    //cout << "drho = " << rho1 << " / " << rho2 << endl;
    }

    if(fabs(Pmax-Pmin)<1e-4)
    {
        rho1 = (rhomax+rhomin)/2;
        rho2 = rho1;
    }


    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;
    //area = maxwell_cons(rho,P,P1);
    //cout << "area newt = " << area << endl;


    rhov[0] = rho1;
    rhov[1] = rho2;
    rhov[2] = P1;
    rhov[3] = P2;
    rhov[4] = u1;
    rhov[5] = u2;
    return rhov;
}

//Function to calculate phase coexistence densities
//Using Newton method
vector<double> dens_newt5(vector<double>& rho, vector<double>& f, vector<double>& P, vector<double>& u, double tol)
{
    std::vector<double> dPdrho(1000), Pf1(1000), Pf2(1000), rhov(6);
    int i, max1, min1, max_b, min_a;
    double rhomax, rhomin, Pmax, Pmin, P1, P2, f1, f2, u1, u2, du1, du2, dP1, dP2, detJ, drho1, drho2;
    double rho1, rho2, area;

    //Calculating vector with first derivatives
    dPdrho = fin_diff_1_vec(rho,P);

    //Find max and min pressure of isotherm inside binodal curve
    max1 = bin_max(dPdrho);
    min1 = bin_min(dPdrho);
    rhomax = rho[max1];
    rhomin = rho[min1];
    Pmax = P[max1];
    Pmin = P[min1];

    max_b = max1-200;
    min_a = min1+200;
    if(max_b-200 < 0) max_b = 0;



    if (Pmin>Pmax)
    {
        min1 = bin_min_seed(dPdrho,max1);
        rhomin = rho[min1];
        Pmin = P[min1];
    }

    if (Pmin<0) Pmin=1e-3;

    for(i=0; i<1000; i++)
    {
        if(Pmin<0) Pf1[i] = P[i] + Pmin;
        else Pf1[i] = P[i] - Pmin;

        Pf2[i] = P[i] - Pmax;
    }

    //cout << "Pmax/min = " << Pmax << " / " << Pmin << endl;
    //cout << "max/min = " << max1 << " / " << min1 << endl;
    //cout << "rho = " << rhomax << " / " << rhomin << endl;
    //Find initial guess for densities
    rho1 = falsi_spline(rho, Pf1, rho[0], rhomax, tol);
    rho2 = falsi_spline(rho, Pf2, rhomin, rho[700], tol);
    //cout << "rho1 = " << rho1 << " / " << rho2 << endl;

    //Solve newton-raphson system
    drho1 = tol+1;
    int counter = 0;
    while(fabs(drho1)>tol || fabs(drho2)>tol)
    {
    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;

    du1 = cspline_deriv1(rho,u,rho1);
    du2 = cspline_deriv1(rho,u,rho2);
    dP1 = cspline_deriv1(rho,P,rho1);
    dP2 = cspline_deriv1(rho,P,rho2);
    detJ = -dP2*du1+dP1*du2;

    drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2);
    drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2);

    rho1 = rho1 + drho1/2;
    rho2 = rho2 + drho2/2;
    //cout << "drho = " << drho1 << " / " << drho2 << endl;
    }

    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;
    //area = maxwell_cons(rho,P,P1);
    //cout << "area newt = " << area << endl;

    rhov[0] = rho1;
    rhov[1] = rho2;
    rhov[2] = P1;
    rhov[3] = P2;
    rhov[4] = u1;
    rhov[5] = u2;
    return rhov;
}


//Function to calculate phase coexistence densities
//Using Newton method
vector<double> dens_newt_seed(vector<double>& rho, vector<double>& f, vector<double>& P, vector<double>& u, double tol, vector<double>& seed)
{
    std::vector<double> dPdrho(1000), Pf1(1000), Pf2(1000), rhov(1000);
    int i, max1, min1;
    double rhomax, rhomin, Pmax, Pmin, P1, P2, f1, f2, u1, u2, du1, du2, dP1, dP2, detJ, drho1, drho2;
    double rho1, rho2, area;

    rho1 = seed[0];
    rho2 = seed[1];

    //Solve newton-raphson system
    drho1 = tol+1;
    int counter = 0;
    while(drho1>tol || drho2>tol)
    {
    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;

    du1 = cspline_deriv1(rho,u,rho1);
    du2 = cspline_deriv1(rho,u,rho2);
    dP1 = cspline_deriv1(rho,P,rho1);
    dP2 = cspline_deriv1(rho,P,rho2);
    detJ = -dP2*du1+dP1*du2;

    drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2);
    drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2);

    rho1 = rho1 + drho1/2;
    rho2 = rho2 + drho2/2;
    //cout << "drho = " << drho1 << " / " << drho2 << endl;
    }

    f1 = cspline(rho,f,rho1);
    f2 = cspline(rho,f,rho2);
    u1 = cspline_deriv1(rho,f,rho1);
    u2 = cspline_deriv1(rho,f,rho2);
    P1 = -f1+rho1*u1;
    P2 = -f2+rho2*u2;
    //area = maxwell_cons(rho,P,P1);
    //cout << "area newt = " << area << endl;

    rhov[0] = rho1;
    rhov[1] = rho2;
    rhov[2] = P1;
    return rhov;
}

vector<double> linear_regression(vector<double>& x, vector<double>& y)
{
    int SIZE, i;
    SIZE = x.size();
    double a, b, sumXY, sumX, sumY, sumX2, sumY2, r, r2;
    std::vector<double> coef(3);
    sumXY = 0;
    sumX = 0;
    sumY = 0;
    sumY2 = 0;
    sumX2 = 0;

    for(i=0; i<SIZE; i++)
    {
        cout << "x / y = " << x[i] << " / " << y[i] << endl;
    }

    for(i=0; i<SIZE; i++)
    {
        sumXY = sumXY + x[i]*y[i];
        sumX  = sumX  + x[i];
        sumY  = sumY  + y[i];
        sumY2 = sumY2 + y[i]*y[i];
        sumX2 = sumX2 + x[i]*x[i];
    cout << "SUM = " << sumXY << " / " << sumX << " / " << sumY << " / " << sumY2 << " / " << sumX2 << endl;
    }

    a = (SIZE*sumXY-sumX*sumY)/(SIZE*sumX2-pow(sumX,2));
    b = (sumY-a*sumX)/SIZE;
    r = (SIZE*sumXY-sumX*sumY)/(pow(SIZE*sumX2-pow(sumX,2),0.5)*pow(SIZE*sumY2-pow(sumY,2),0.5));
    r2 = pow(r,2);

    cout << "a / b = " << a << " / " << b << endl;
    cout << "R2 = " << r2 << endl;
    coef[0] = a;
    coef[1] = b;
    coef[2] = r2;

    return coef;
}

double linear_regression_angular(vector<double>& x, vector<double>& y)
{
    int SIZE, i;
    SIZE = x.size();
    double a, b, sumXY, sumX, sumY, sumX2, sumY2, r, r2;
    std::vector<double> coef(3);
    sumXY = 0;
    sumX = 0;
    sumY = 0;
    sumY2 = 0;
    sumX2 = 0;

    /*
    for(i=0; i<SIZE; i++)
    {
        cout << "x / y = " << x[i] << " / " << y[i] << endl;
    }
    */

    for(i=0; i<SIZE; i++)
    {
        sumXY = sumXY + x[i]*y[i];
        sumX  = sumX  + x[i];
        sumY  = sumY  + y[i];
        sumY2 = sumY2 + y[i]*y[i];
        sumX2 = sumX2 + x[i]*x[i];
    //cout << "SUM = " << sumXY << " / " << sumX << " / " << sumY << " / " << sumY2 << " / " << sumX2 << endl;
    }

    a = (SIZE*sumXY-sumX*sumY)/(SIZE*sumX2-pow(sumX,2));
    b = (sumY-a*sumX)/SIZE;
    r = (SIZE*sumXY-sumX*sumY)/(pow(SIZE*sumX2-pow(sumX,2),0.5)*pow(SIZE*sumY2-pow(sumY,2),0.5));
    r2 = pow(r,2);

    //cout << "a / b = " << a << " / " << b << endl;
    //cout << "R2 = " << r2 << endl;

    return a;
}

#endif // NUMERICAL_H_INCLUDED
