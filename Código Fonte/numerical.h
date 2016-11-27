#ifndef NUMERICAL_H_INCLUDED
#define NUMERICAL_H_INCLUDED
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "interpolation_util.h"
#include <math.h>

//Trapezoidal method to calculate area of integral, , uses discrete data points
//Outputs area under function and given interval
double trapezoidal_rule(vector<double>& x, vector<double>& y, int SIZE)
{
int i = 0;
//int SIZE = x.size();
double sum = 0;
double h = (x[SIZE]-x[0])/SIZE;
//cout << "y = " << y[0] << " / " << 0 << endl;
sum = sum + y[0];
for(i=1; i<SIZE; i++)
{
sum = sum + 2*y[i];
//cout << "y = " << y[i] << " / " << i << " / sum parcial = " << sum << endl;
}
//cout << "y = " << y[SIZE] << " / " << SIZE << endl;
sum = sum + y[SIZE];
//cout << "sum = " << sum << " | h/2 = " << h/2;
sum = sum * h/2;
return sum;
}


//Simpson method to calculate area of integral, uses discrete data points
//Outputs area under function and given interval
double simpson_rule(vector<double>& x, vector<double>& y, int SIZE)
{
int i = 0;
//int SIZE = x.size();
double sum = 0;
double h = (x[SIZE]-x[0])/SIZE;
for(i=1; i<SIZE; i++)
{
if(i%2 > 0) sum = sum + 4*y[i];
else sum = sum + 2*y[i];
}
sum = sum + y[0] + y[SIZE];
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
trap_area = simpson_rule(x_subv,y_subv, b-a);
Area = trap_area - (y[b]+y[a])/2*(x[b]-x[a]);
return Area;
}


//Finite difference method to calculate approximate first derivative
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
y1[SIZE] = (y[SIZE]-y[SIZE2])/(x[SIZE]-x[SIZE2]);

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
std::vector<double> y1(SIZE); //Vector to output first derivatives
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
double phase_coexistence_e(double T, VectorXd P, VectorXd u, VectorXd rho, int phase)
{
    int SIZE, i, SIZE2;
    double rho_liq, rho_vap, rho1, drhol, drhov, dPl, dPv, dul, duv, Pliq, Pvap, uliq, uvap, detJ, cond_rhov, cond_rhol, tol;
    SIZE = P.size(); //Measure size of vectors
    SIZE2 = SIZE-1;

    //std::vector<double> *rho2, *dP_drho2, *du_drho2, *P2, *u2;

    VectorXd dP_drho(SIZE), du_drho(SIZE);
    VectorXd F(2);
    MatrixXd J(2,2);

    cout << "before cond " << endl;
    tol = 1e-3;
    cond_rhol = tol+1;
    cond_rhov = tol+1;

    //Initial guess for solution
    cout << "before guess" << endl;
    Pliq = P(990);
    Pvap = P(10);
    uliq = u(990);
    uvap = u(10);
    rho_liq = rho(990);
    rho_vap = rho(10);

    //Pressure and chemical potential derivatives for density
    cout << "before deriv" << endl;
    dP_drho(0) = (P(1)-P(0))/(rho(1)-rho(0));
    for(i=1; i<SIZE-1; i++)
    {
    dP_drho(i) = (P(i+1)-P(i-1))/(rho(i+1)-rho(i-1));
    }
    dP_drho(SIZE) = (P(SIZE)-P(SIZE2))/(rho(SIZE)-rho(SIZE2));

    du_drho(0) = (u(1)-u(0))/(rho(1)-rho(0));
    for(i=1; i<SIZE-1; i++)
    {
    du_drho(i) = (u(i+1)-u(i-1))/(rho(i+1)-rho(i-1));
    }
    du_drho(SIZE) = (u(SIZE)-u(SIZE2))/(rho(SIZE)-rho(SIZE2));

    //dP_drho = fin_diff_1_vec(*rho, *P);
    //du_drho = fin_diff_1_vec(*rho, *u);
    cout << "before iter" << endl;

    std::vector<double> dP_drho2(1000), du_drho2(1000), rho2(1000), P2(1000), u2(1000);

for(i=0;i<1001;i++)
{
    dP_drho2[i] = dP_drho(i);
    du_drho2[i] = du_drho(i);
    rho2[i] = rho(i);
    P2[i] = P(i);
    u2[i] = u(i);
}

    cout << "dPdrho = " << dP_drho2[100] << " " << dP_drho[100] << endl;
/*
        vector<double> dP_drho2(dP_drho.data(), dP_drho.data() + dP_drho.rows() * dP_drho.cols());
        vector<double> du_drho2(du_drho.data(), du_drho.data() + du_drho.rows() * du_drho.cols());
        vector<double> rho2(rho.data(), rho.data() + rho.rows() * rho.cols());
        vector<double> P2(P.data(), P.data() + P.rows() * P.cols());
        vector<double> u2(u.data(), u.data() + u.rows() * u.cols());
*/

cout << "before conditions" << endl;
    while(cond_rhol > tol || cond_rhov > tol)
    {
    //Objective function
    F(0) = Pvap - Pliq;
    F(1) = uvap - uliq;

    cout << "before spline" << endl;
    dPl = cspline(rho2,dP_drho2,rho_liq);
    dPv = cspline(rho2,dP_drho2,rho_vap);
    dul = cspline(rho2,du_drho2,rho_liq);
    duv = cspline(rho2,du_drho2,rho_vap);


    cout << "before J" << endl;
    //jacobian Matrix
    J(0,0) = -dPl;
    J(0,1) = dPv;
    J(1,0) = -dul;
    J(1,1) = duv;
    detJ = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    cout << "J = " << J << endl;
    cout << "detJ = " << detJ << endl;

    //Solving J^-1 * F
    drhol = -duv/detJ*F(0)+dPv/detJ*F(1);
    drhov = -dul/detJ*F(0)+dPl/detJ*F(1);

    rho_liq = rho_liq + drhol;
    rho_vap = rho_vap + drhov;

    Pliq = cspline(rho2,P2,rho_liq);
    Pvap = cspline(rho2,P2,rho_vap);
    uliq = cspline(rho2,u2,rho_liq);
    uvap = cspline(rho2,u2,rho_vap);

    cond_rhol = fabs(drhol/rho_liq);
    cond_rhov = fabs(drhov/rho_vap);
    cout << "cond_rhol = " << cond_rhol << " " << cond_rhov << endl;
    cout << "liq / vap = " << rho_liq << " / " << rho_vap << endl;
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



#endif // NUMERICAL_H_INCLUDED
