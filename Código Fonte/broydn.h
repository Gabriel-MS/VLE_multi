#ifndef BROYDN_H_INCLUDED
#define BROYDN_H_INCLUDED

#include <math.h>
#include "nrutil.h"
#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6

#define ALF 1.0e-4

#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);\
free_matrix(r,1,n,1,n);free_matrix(qt,1,n,1,n);free_vector(p,1,n);\
free_vector(g,1,n);free_vector(fvcold,1,n);free_vector(d,1,n);\
free_vector(c,1,n);return;}


void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
double *f, double stpmax, int *check, double (*func)(double []))
{
int i;
double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
test,tmplam;
*check=0;
for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
sum=sqrt(sum);
if (sum > stpmax)
for (i=1;i<=n;i++) p[i] *= stpmax/sum;
for (slope=0.0,i=1;i<=n;i++)
slope += g[i]*p[i];
if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
test=0.0;
for (i=1;i<=n;i++) {
temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
if (temp > test) test=temp;
}
alamin=TOLX/test;
alam=1.0;
for (;;) {
for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
*f=(*func)(x);
if (alam < alamin) {
for (i=1;i<=n;i++) x[i]=xold[i];
*check=1;
return;
} else if (*f <= fold+ALF*alam*slope) return;
else { Backtrack.
if (alam == 1.0)
tmplam = -slope/(2.0*(*f-fold-slope));
else {
rhs1 = *f-fold-alam*slope;
rhs2=f2-fold-alam2*slope;
a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
if (a == 0.0) tmplam = -slope/(2.0*b);
else {
disc=b*b-3.0*a*slope;
if (disc < 0.0) tmplam=0.5*alam;
else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
else tmplam=-slope/(b+sqrt(disc));
}
if (tmplam > 0.5*alam)
tmplam=0.5*alam; λ ≤ 0.5λ1.
}
}
alam2=alam;
f2 = *f;
alam=FMAX(tmplam,0.1*alam); λ ≥ 0.1λ1.
}
}

void fdjac(int n, double x[], double fvec[], double **df,
void (*vecfunc)(int, double [], double []))
{
int i,j;
double h,temp,*f;
f=vector(1,n);
for (j=1;j<=n;j++) {
temp=x[j];
h=EPS*fabs(temp);
if (h == 0.0) h=EPS;
x[j]=temp+h;
h=x[j]-temp;
(*vecfunc)(n,x,f);
x[j]=temp;
for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
}
free_vector(f,1,n);
}

extern int nn;
extern double *fvec;
extern void (*nrfuncv)(int n, double v[], double f[]);
double fmin(double x[])
{
int i;
double sum;
(*nrfuncv)(nn,x,fvec);
for (sum=0.0,i=1;i<=nn;i++) sum += SQR(fvec[i]);
return 0.5*sum;
}

void qrdcmp(double **a, int n, double *c, double *d, int *sing)
{
int i,j,k;
double scale,sigma,sum,tau;
*sing=0;
for (k=1;k<n;k++) {
scale=0.0;
for (i=k;i<=n;i++) scale=FMAX(scale,fabs(a[i][k]));
if (scale == 0.0) {
*sing=1;
c[k]=d[k]=0.0;
} else {
for (i=k;i<=n;i++) a[i][k] /= scale;
for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
sigma=SIGN(sqrt(sum),a[k][k]);
a[k][k] += sigma;
c[k]=sigma*a[k][k];
d[k] = -scale*sigma;
for (j=k+1;j<=n;j++) {
for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
tau=sum/c[k];
for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
}
}
}
d[n]=a[n][n];
if (d[n] == 0.0) *sing=1;
}

void rsolv(double **a, int n, double d[], double b[])
{
int i,j;
double sum;
b[n] /= d[n];
for (i=n-1;i>=1;i--) {
for (sum=0.0,j=i+1;j<=n;j++) sum += a[i][j]*b[j];
b[i]=(b[i]-sum)/d[i];
}
}

void qrsolv(double **a, int n, double c[], double d[], double b[])
{
void rsolv(double **a, int n, double d[], double b[]);
int i,j;
double sum,tau;
for (j=1;j<n;j++) {
for (sum=0.0,i=j;i<=n;i++) sum += a[i][j]*b[i];
tau=sum/c[j];
for (i=j;i<=n;i++) b[i] -= tau*a[i][j];
}
rsolv(a,n,d,b);
}

void rotate(double **r, double **qt, int n, int i, double a, double b)
{
int j;
double c,fact,s,w,y;
if (a == 0.0) {
c=0.0;
s=(b >= 0.0 ? 1.0 : -1.0);
} else if (fabs(a) > fabs(b)) {
fact=b/a;
c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
s=fact*c;
} else {
fact=a/b;
s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
c=fact*s;
}
for (j=i;j<=n;j++) {
y=r[i][j];
w=r[i+1][j];
r[i][j]=c*y-s*w;
r[i+1][j]=s*y+c*w;
}
for (j=1;j<=n;j++) {
y=qt[i][j];
w=qt[i+1][j];
qt[i][j]=c*y-s*w;
qt[i+1][j]=s*y+c*w;
}
}

void qrupdt(double **r, double **qt, int n, double u[], double v[])
{
void rotate(double **r, double **qt, int n, int i, double a, double b);
int i,j,k;
for (k=n;k>=1;k--) {
if (u[k]) break;
}
if (k < 1) k=1;
for (i=k-1;i>=1;i--) {
rotate(r,qt,n,i,u[i],-u[i+1]);
if (u[i] == 0.0) u[i]=fabs(u[i+1]);
else if (fabs(u[i]) > fabs(u[i+1]))
u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
else u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
}
for (j=1;j<=n;j++) r[1][j] += u[1]*v[j];
for (i=1;i<k;i++) rotate(r,qt,n,i,r[i][i],-r[i+1][i]);
}

/*
void vecfunc(int n, double x[], double fvec[])
{
    fvec = dens_solver(double x)
}
*/

int nn;
double *fvec;
void (*nrfuncv)(int n, double v[], double f[]);

//Broyden's method to solve nonlinear systems of equation
void broydn(double x[], int n, int *check,
void (*vecfunc)(int, double [], double []))
{

void fdjac(int n, double x[], double fvec[], double **df,
void (*vecfunc)(int, double [], double []));
double fmin(double x[]);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
double *f, double stpmax, int *check, double (*func)(double []));
void qrdcmp(double **a, int n, double *c, double *d, int *sing);
void qrupdt(double **r, double **qt, int n, double u[], double v[]);
void rsolv(double **a, int n, double d[], double b[]);
int i,its,j,k,restrt,sing,skip;
double den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
double *g,*p,**qt,**r,*s,*t,*w,*xold;
c=vector(1,n);
d=vector(1,n);
fvcold=vector(1,n);
g=vector(1,n);
p=vector(1,n);
qt=matrix(1,n,1,n);
r=matrix(1,n,1,n);
s=vector(1,n);
t=vector(1,n);
w=vector(1,n);
xold=vector(1,n);
fvec=vector(1,n);
nn=n;
nrfuncv=vecfunc;
f=fmin(x);
test=0.0;
for (i=1;i<=n;i++)
if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
if (test < 0.01*TOLF) {
*check=0;
FREERETURN
}

for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
stpmax=STPMX*FMAX(sqrt(sum),(double)n);
restrt=1;
for (its=1;its<=MAXITS;its++) {
if (restrt) {
fdjac(n,x,fvec,r,vecfunc);
qrdcmp(r,n,c,d,&sing);
if (sing) nrerror("singular Jacobian in broydn");
for (i=1;i<=n;i++) {
for (j=1;j<=n;j++) qt[i][j]=0.0;
qt[i][i]=1.0;
}
for (k=1;k<n;k++) {
if (c[k]) {
for (j=1;j<=n;j++) {
sum=0.0;
for (i=k;i<=n;i++)
sum += r[i][k]*qt[i][j];
sum /= c[k];
for (i=k;i<=n;i++)
qt[i][j] -= sum*r[i][k];
}
}
}
for (i=1;i<=n;i++) {
r[i][i]=d[i];
for (j=1;j<i;j++) r[i][j]=0.0;
}
} else {
for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
for (i=1;i<=n;i++) {
for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
t[i]=sum;
}
skip=1;
for (i=1;i<=n;i++) {
for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
w[i]=fvec[i]-fvcold[i]-sum;
if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;

else w[i]=0.0;
}
if (!skip) {
for (i=1;i<=n;i++) {
for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
t[i]=sum;
}

for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
for (i=1;i<=n;i++) s[i] /= den;
qrupdt(r,qt,n,t,s);
for (i=1;i<=n;i++) {
if (r[i][i] == 0.0) nrerror("r singular in broydn");
d[i]=r[i][i];
}
}
}

 for (i=1;i<=n;i++) {
for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
g[i]=sum;
}
for (i=n;i>=1;i--) {
for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
g[i]=sum;
}
for (i=1;i<=n;i++) {
xold[i]=x[i];
fvcold[i]=fvec[i];
}
fold=f; Store f.
for (i=1;i<=n;i++) {
for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
p[i] = -sum;
}
rsolv(r,n,d,p);
lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
test=0.0;
for (i=1;i<=n;i++)
if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
if (test < TOLF) {
*check=0;
FREERETURN
}
if (*check) {
if (restrt) FREERETURN Failure;
{
test=0.0;
den=FMAX(f,0.5*n);
for (i=1;i<=n;i++) {
temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
if (temp > test) test=temp;
}
if (test < TOLMIN) FREERETURN
else restrt=1;
}
} else {
restrt=0;
test=0.0;
for (i=1;i<=n;i++) {
temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
if (temp > test) test=temp;
}
if (test < TOLX) FREERETURN
}
}
nrerror("MAXITS exceeded in broydn");
FREERETURN
}





#endif // BROYDN_H_INCLUDED
