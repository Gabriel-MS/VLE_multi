#ifndef NUMERICAL_H_INCLUDED
#define NUMERICAL_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
//#include "nrutil.h" //using numerical recipes c++
#include "interpolation_util.h"
#include <math.h>

//Finite difference method to calculate approximate first derivative
//x and y vectors must have the same size
//Outputs first derivatives in vector
vector<double> fin_diff_1_vec(const vector<double>& x, const vector<double>& y)
{
  int SIZE, i;

  SIZE = x.size();

  std::vector<double> y1(SIZE); //Vector to output first derivatives

  y1[0] = y[0]/x[0]; //Wrong, fix

  for(i=1; i<SIZE; i++)
  {
    y1[i] = (y[i]-y[i-1])/(x[i]-x[i-1]);
  }

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


//Cubic spline, outputs 2nd derivative to use in interpolation function
//end1 and endn are the specified first derivatives boundaries
//at the begin and end respectively
/*                     UNDER DEVELOPMENT
double cubic_spline(double x[], double y[], int n, double end1, double endn, double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;

  u=vector(1,n-1);
  if (yp1 > 0.99e30) y2[1]=u[1]=0.0;
    else
    {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }


  for (i=2;i<=n-1;i++)
  {
  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p=sig*y2[i-1]+2.0;
  y2[i]=(sig-1.0)/p;
  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  if (ypn > 0.99e30) qn=un=0.0;
    else
    {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }

  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

  for (k=n-1;k>=1;k--)  y2[k]=y2[k]*y2[k+1]+u[k];

  free_vector(u,1,n-1);
}
*/


//Cubic spline interpolation, interpolates single values in cubic splines
double cspline(const vector<double>& X, const vector<double>& Y, double x)
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
vector<double> cspline_vec(const vector<double>& X, const vector<double>& Y, vector<double>& x)
{
  int SIZE, i;

  tk::spline s;
  s.set_points(X,Y);

  SIZE = X.size();

  std::vector<double> y1(SIZE); //Vector to output first derivatives

  for(i=1; i<SIZE; i++)
  {
    y1[i] = s(x[i]);
  }

return y1;
}


//Cubic spline derivative, outputs values of 1st derivatives in cubic splines
double cspline_deriv1(const vector<double>& X, const vector<double>& Y, double x)
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

  SIZE = X.size();

  std::vector<double> y1(SIZE); //Vector to output first derivatives

  for(i=1; i<SIZE; i++)
  {
    y1[i] = s.deriv(1,x[i]);
  }

return y1;
}

#endif // NUMERICAL_H_INCLUDED
