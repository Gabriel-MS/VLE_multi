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
double trapezoidal_interval(vector<double>& x, vector<double>& y, int SIZE, int a, int b)
{
int i = 0;
double sum = 0;
double h = (x[SIZE-1]-x[0])/SIZE;
sum = sum + y[a];
for(i=a; i<b+1; i++)
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
double simpson_interval(vector<double>& x, vector<double>& y, int SIZE, int a, int b)
{
int i = 0;
//int SIZE = x.size();
double sum = 0;
double h = (x[SIZE-1]-x[0])/SIZE;
for(i=a; i<b+1; i++)
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
y2[0] = (y[1]- 2*y[0])/(pow(x[1]-x[0],2)); //Wrong, fix
for(i=1; i<SIZE; i++)
{
y2[i] = (y[i+1]-2*y[i]+y[i-1])/(pow(x[i]-x[i-1],2));
}
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
int i=100;
do
{
 pmax_cond = dP_dV[i];
 i++;
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
//cout << "j = " << j << " / pmin = " << pmin_cond << " / rho = " << rho_vec_out[j] << endl;
}while(pmin_cond>0);
    return j+1;
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
    cout << "falsi a b c = " << a << " / " << b << " / " << c << endl;
    if(ya*yc<0) b = c;
    else a = c;
    //cout << "AFTER falsi a b c = " << a << " / " << b << " / " << c << endl;
    //cin >> i;
    }

    return c;
}

//Calculates Area of given isotherm (helmholtz energy density f, density rho and chemical potential u)
//All inputs should be provided dimensionless
//Uses Cubic Spline to interpolate data points
double Maxwell_cons(vector<double>& rho, vector<double>& P, double P1)
{
    cout << "Maxwell in" << endl;
    std::vector<double> dPdV(1000), rho_a1(100), rho_a2(100), P_a1(100), P_a2(100), Pf(1000);
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
    //cout << "rho_max = " << rho_max << " / " << rho_min << endl;
    P_max = P[max1];
    P_min = P[min1];
    //cout << "P1 = " << P1 << endl;

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

    //cout << "Pmin = " << P[min1] << endl;
    //cout << "Pmax = " << P[max1] << endl;

    //Find density roots for pressure using Regula Falsi method
    rho1 = falsi_spline(rho, Pf, rho[0], rho_max, 1e-3);
    cout << "rho1 = " << rho1 << endl;
    rho2 = falsi_spline(rho, Pf, rho_max, rho_min, 1e-3);
    //cout << "rho2 = " << rho2 << endl;
    rho3 = falsi_spline(rho, Pf, rho_min, rho[700], 1e-3);
    //cout << "rho3 = " << rho3 << endl;

    //Find new values of pressure and density between rho1 and rho3
    rhoz1 = rho1;
    rhoz2 = rho2;
    h1 = (rho2-rho1)/100;
    h2 = (rho3-rho2)/100;

    rho_a1[0] = rho1;
    rho_a1[99] = rho2;
    rho_a2[0] = rho2;
    rho_a2[99] = rho3;

    P_a1[0] = P1;
    P_a1[99] = P1;
    P_a2[0] = P1;
    P_a2[99] = P1;

    for(i=1; i<99; i++)
    {
    rho_a1[i] = rhoz1;
    rho_a2[i] = rhoz2;
    rhoz1 = rhoz1 + h1;
    rhoz2 = rhoz2 + h2;
    }

    for(i=1; i<99; i++)
    {
    P_a1[i] = cspline(rho, P, rho_a1[i]);
    P_a2[i] = cspline(rho, P, rho_a2[i]);
    //cout << "a1 = " << rho_a1[i] << " / " << P_a1[i] << endl;
    }

    //Calculates Maxwell area
    Area1 = simpson_rule(rho_a1, P_a1);
    Area1 = Area1 - P1*(rho2-rho1);

    Area2 = P1*(rho3-rho2);
    Area2 = Area2 - simpson_rule(rho_a2, P_a2);

    Area = Area1 - Area2;
    Area = fabs(Area);

    //cout << "P1 = " << P1 << endl;
    //cout << "rho1 = " << rho1 << endl;
    //cout << "rho2 = " << rho2 << endl;
    //cout << "rho3 = " << rho3 << endl;
    //cout << "Area1 pre = " << simpson_rule(rho_a1, P_a1) << endl;
    //cout << "Areas = " << Area1 << " / " << Area2 << endl;
    //cout << "Area = " << Area << endl;
    cout << "Maxwell out" << endl;
    //cin >> i;
    return Area;
}

//Function to calculate density of coexistent phases using maxwell construction
//Input is dimensionless isotherm
vector<double> dens_maxwell(vector<double>& rho, vector<double>& P)
{
    int i, SIZE;
    std::vector<double> dens(3);
    std::vector<double> dPdV(1000), Pf(1000);
    double rho1, rho2, rho_max, rho_min, P_max, P_min, P0, P1, P2, Pz;
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

    cout << "before brent" << endl;

    //Brent method to use Maxwell Construction Method
    Pz = zbrentm(Maxwell_cons, P1*0.95, P2*1.05, 1e-5, rho, P);
    //cout << "Pz = " << Pz << endl;
    cout << "after brent" << endl;

    //Adjust isotherm to find roots to pressure
    for(i=0; i<SIZE+1; i++)
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
    cout << "rho1 = " << rho1;
    cout << "rho2 = " << rho2;
    return dens;
}

#endif // NUMERICAL_H_INCLUDED
