#ifndef BICUBIC_H_INCLUDED
#define BICUBIC_H_INCLUDED

#include <vector>

using namespace std;

void spline(double* x, double* y, double yp1, double ypn, int n, double* y2)
{
    //int n = x.size();
    //std::vector<double> u(n-1), y2(n);
    std::vector<double> u(n-1);
    n = n - 1;

int i,k;
double p,qn,sig,un;

if (yp1 > 0.99e30)
y2[0]=u[0]=0.0;
else {
y2[0] = -0.5;
u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
}


for (i=2;i<=n-1;i++) {
sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
p=sig*y2[i-1]+2.0;
y2[i]=(sig-1.0)/p;
u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
}
if (ypn > 0.99e30)
qn=un=0.0;
else {
qn=0.5;
un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
}

y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
for (k=n-1;k>=1;k--)
y2[k]=y2[k]*y2[k+1]+u[k];

//return y2;

}

void splint(double* xa, double* ya, double* y2a, double x, int n, double *y)
{
    n = n-1;


int klo,khi,k;
double h,b,a;
klo=0;
khi=n;
while (khi-klo > 1) {
k=(khi+klo) >> 1;
if (xa[k] > x) khi=k;
else klo=k;
}
h=xa[khi]-xa[klo];

a=(xa[khi]-x)/h;
b=(x-xa[klo])/h;
*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void splie2(double* x1a, double* x2a, int m, int n, double **ya, double **y2a)
{
  m = m - 1;
  n = n - 1;
  int j;

  for (j=0;j<=m;j++)
  spline(x2a,ya[j],1.0e30,1.0e30,n,y2a[j]);
}

void splin2(double* x1a, double* x2a, double **ya, double **y2a, int m, int n,
double x1, double x2, double *y)
{
  int j;
  double ytmp[m];
  double yytmp[m];
  double yytmpval;


  for (j=0;j<=m-1;j++){
  splint(x2a,ya[j],y2a[j],x2,n,&yytmpval);
  yytmp[j] = yytmpval;
  }
  spline(x1a,yytmp,1.0e30,1.0e30,m,ytmp);
  splint(x1a,yytmp,ytmp,x1,m,y);
}

double bicubic_s_int(vector<double>& x1, vector<double>& x2, double x1int, double x2int, double **ym)
{
    //x1 makes reference to rows, size n1
    //x2 makes reference to cols, size n2
    double y; //y value to be calculated by interpolation
    int n1 = x1.size();
    int n2 = x2.size();
    double x1a[n1], x2a[n2];

    for(int i=0;i<n1;i++)
    {
        x1a[i] = x1[i];
    }

    for(int i=0;i<n2;i++)
    {
        x2a[i] = x2[i];
    }

    double** y2; //Second derivative 2d array
    y2 = new double *[n1];
    for(int k = 0; k <n1; k++)
        y2[k] = new double[n2];

    splie2(x1a,x2a,n1,n2,ym,y2); //Generates 2nd derivative array y
    splin2(x1a,x2a,ym,y2,n1,n2,x1int,x2int,&y);

    for(int i=0; i<n1; i++)
    {
        delete[] y2[i];
    }
    delete[] y2;

    return y;
}






#endif // BICUBIC_H_INCLUDED
