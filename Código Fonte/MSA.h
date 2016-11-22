#ifndef MSA_H_INCLUDED
#define MSA_H_INCLUDED

#include <cmath>
#include <complex>
#include <vector>
/*
complex operator + (complex a, int b)
{
    return complex(a.real() * b, a.imag() * b);
}
*/

double script_f(double tstar, double y, double xlam)
{
      std::complex<double> l, fj, f1j, m, m1, s, s1, s2, fi, f1i, ej, e1j, e2j, n, n1, a, a1, sum1, sum2, cj, c1;
      std::complex<double> k0, k1, k2, c_0, c_1, c_2, c2, n2, d2, d, d1;
      double f, arg, yp, ym, u1, aa1, zcomp1, aa2, u2, zcomp2, zcomp, fac, dfac, pi, pyf, s3;
      vector<double> kr(3), ki(3);
      int i, j;
      pi = 3.14159265359796;
      std::complex<double> c(cos(2*pi/3),sin(2*pi/3));
      std::complex<double> sum(0,0);

      f = 3+3*y-y*y;
      arg=pow((1+2*pow(y,4)/(f*f)),0.5);
      yp= pow((arg+1),(1/3));
      ym=-pow((arg-1),(1/3));

        //for(j=0;j<3;j++)
        //{
         //cj = pow(c,j);
         pyf = pow((2*y*f),(1/3));
         c_0 = pow(c,0);
         c_1 = pow(c,1);
         c_2 = pow(c,2);
         k0 = -2*y+pyf*(yp*c_0+ym/c_0)/(1-y);
         k1 = -2*y+pyf*(yp*c_1+ym/c_1)/(1-y);
         k2 = -2*y+pyf*(yp*c_2+ym/c_2)/(1-y);

         kr[0] = k0.real();
         kr[1] = k1.real();
         kr[2] = k2.real();
         ki[0] = k0.imag();
         ki[1] = k1.imag();
         ki[2] = k2.imag();


         vector< complex<double> > t, tm;
         for(i=0; i<3; i++)
         {
         t.push_back( complex<double>( kr[i], ki[i] ) );
         tm.push_back( complex<double>( -kr[i], -ki[i] ) );
         }

         cout << "\n===================" << endl;
         cout << "t = " << tstar << endl;
         cout << "y = " << y << endl;
         cout << "t0 = " << t[0] << endl;
         cout << "t1 = " << t[1] << endl;
         cout << "t2 = " << t[2] << endl;
         cin >> i;

         //std::complex<double>   t[2] = {-2*y+pow((2*y*f),(1/3))*(yp*cj+ym/cj)/(1-y)};
        //}
//=================first order contribution===========================
      for(i=0;i<3;i++)
        {
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         l = (1.0+0.5*y)*t[i]+1.0+2.0*y;
         sum = sum + l*(xlam*t[i]-1.0)*exp((xlam-1.0)*t[i])/(t[i]*s1);
        }

      u1= -1 - 12*y*sum.real();
      aa1=u1;

      sum = sum - sum;
      for(i=0;i<3;i++)
        {
         s1 = 3.0*pow((1-y),2.0)*t[i]*t[i]+12.0*y*(1.0-y)*t[i]+18.0*y*y;
         s2 = 6.0*pow((1.0-y),2.0)*t[i]+12.0*y*(1.0-y);
         a = xlam*(1.0-y)*(1.0+2.0*y)*t[i]*t[i]*t[i]*t[i] +(xlam*pow((1.0+2.0*y),2.0)-(1.0-y)*(1.0+2.0*y))*t[i]*t[i]*t[i] -pow((1.0+2.0*y),2.0)*t[i]*t[i];
         a1 = 4.0*xlam*(1.0-y)*(1.0+2.0*y)*t[i]*t[i]*t[i] +3.0*(xlam*pow((1.0+2.0*y),2.0)-(1.0-y)*(1.0+2.0*y))*t[i]*t[i] -2.0*pow((1.0+2.0*y),2.0)*t[i];
         sum = sum + ((a1+(xlam-1.0)*a)/(s1*s1)-s2*a/(s1*s1*s1))*exp(t[i]*(xlam-1.0));
        }

      zcomp1= -12*y*sum.real();

//==================second order contribution=================

      sum1 = sum1 - sum1;
      sum2 = sum2 - sum2;

      for(i=0;i<3;i++)
        {
         s = pow((1-y),2)*t[i]*t[i]*t[i]+6*y*(1-y)*t[i]*t[i] + 18*y*y*t[i]-12*y*(1+2*y);
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         s2 = 6*pow((1-y),2)*t[i]+12*y*(1-y);
         m = xlam*pow(t[i],5)-pow(t[i],4);
         c1 = 5.0*xlam*pow(t[i],4.0)-4.0*pow(t[i],3.0);
         m1 = c1+(xlam-1)*m;
         fi = m/(s1*s1);
         f1i = m1/(s1*s1)-m*s2/(s1*s1*s1);

         sum1 = sum1+(m*f1i/(s*s)-(m1*s-2.0*m*s1)*fi/(s*s*s));

         for(j=0;j<3;j++)
         {
            s1 = 3*pow((1-y),2)*t[j]*t[j]+12*y*(1-y)*t[j]+18*y*y;
            s2 = 6*pow((1-y),2)*t[j]+12*y*(1-y);
            m = xlam*pow(t[j],5)-pow(t[j],4);
            c1 = 50*xlam*pow(t[j],4.0)-4.0*pow(t[j],3.0);
            m1 = c1+(xlam-1)*m;
            fj = m/(s1*s1);
            f1j = m1/(s1*s1)-m*s2/(s1*s1*s1);

            sum2 = sum2+(f1i*f1j-(fi*f1j+f1i*fj)/(t[i]+t[j]) +2.0*fi*fj/pow((t[i]+t[j]),2.0)) *exp((t[i]+t[j])*(xlam-1))/(t[i]+t[j]);
         }
        }

      u2 = -12*y*pow((1-y),8)*(sum1+sum2).real();
      aa2 = 0.5*u2;
//====================================================
      sum1 = sum1 - sum1;
      sum2 = sum2 - sum2;

      for(i=0;i<3;i++)
        {
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         s2 = 6*pow((1-y),2)*t[i]+12*y*(1-y);
         m = xlam*pow(t[i],5)-pow(t[i],4);
         c1 = 5.0*xlam*pow(t[i],4.0)-4.0*pow(t[i],3.0);
         m1 = c1+(xlam-1)*m;
         fi = m/(s1*s1);
         f1i = m1/(s1*s1)-m*s2/(s1*s1*s1);

         s = pow((1-y),2)*tm[i]*tm[i]*tm[i]+6*y*(1-y)*tm[i]*tm[i] + 18*y*y*tm[i]-12*y*(1+2*y);
         s1 = 3*pow((1-y),2)*tm[i]*tm[i]+12*y*(1-y)*tm[i]+18*y*y;
         c = xlam*pow(tm[i],5)-pow(tm[i],4);
         d = pow((1.0-y),3.0)*tm[i]*tm[i]*tm[i] -6.0*y*(1.0-y)*(3.0+y)*tm[i]*tm[i] -18.0*y*y*(7.0+y)*tm[i] +12.0*y*(3.0+19.0*y+2.0*y*y);
         n = c*d;
         c1 = 5.0*xlam*pow(tm[i],4.0)-4.0*pow(tm[i],3.0);
         d1 = 3.0*pow((1-y),3.0)*tm[i]*tm[i] -12.0*y*(1.0-y)*(3.0+y)*tm[i] -18.0*y*y*(7.0+y);
         n1 = c1*d+c*d1+(xlam-1)*c*d;

         sum1 = sum1+(n*f1i/(s*s*s)-(n1*s-3.0*n*s1)*fi/(s*s*s*s));

         for(j=0;j<3;j++)
            {
            c = xlam*pow(t[j],5.0)-pow(t[j],4.0);
            d = pow((1-y),3)*t[j]*t[j]*t[j] -6*y*(1-y)*(3+y)*t[j]*t[j] -18*y*y*(7+y)*t[j] +12*y*(3+19*y+2*y*y);
            n = c*d;
            s1 = 3*pow((1-y),2)*t[j]*t[j]+12*y*(1-y)*t[j]+18*y*y;
            ej = n/(s1*s1*s1);

            s2 = 6*pow((1-y),2)*t[j]+12*y*(1-y);
            c1 = 5.0*xlam*pow(t[j],4.0)-4.0*pow(t[j],3.0);
            d1 = 3*pow((1-y),3)*t[j]*t[j] -12*y*(1-y)*(3+y)*t[j] -18*y*y*(7+y);
            n1 = c1*d+c*d1+(xlam-1)*c*d;
            e1j = n1/(s1*s1*s1)-3.0*n*s2/(2.0*pow(s1,4.0));

            c2=20.0*xlam*pow(t[j],3.0)-12.0*pow(t[j],2.0);
            d2=6.0*pow((1-y),3)*t[j] -12*y*(1-y)*(3+y);
            n2 = c2*d+c*d2+pow((xlam-1),2.0)*c*d+2.0*c1*d1 +2.0*(xlam-1.0)*c1*d+2.0*(xlam-1.0)*c*d1;
            s3 = 6.0*pow((1-y),2);
            e2j = n2/(s1*s1*s1)-3.0*n1*s2/pow(s1,4.0) +n*(3.0*s2*s2-s1*s3)/pow(s1,5.0);

            sum2=sum2+(0.5*f1i*e2j -(fi*e2j+2.0*f1i*e1j)/(2.0*(t[i]+t[j])) +(2.0*fi*e1j+f1i*ej)/pow((t[i]+t[j]),2.0) -(3.0*fi*ej)/pow((t[i]+t[j]),3.0)) *exp((t[i]+t[j])*(xlam-1))/(t[i]+t[j]);
            }
        }

      zcomp2 = -6*y*pow((1-y),7)*(sum1+sum2).real();

      zcomp = zcomp1 + zcomp2/tstar;

      fac=(aa1+aa2/tstar)/y;

      dfac=zcomp/(y*y)-fac/y;

      return fac;
}

double script_df(double tstar, double y, double xlam)
{
      std::complex<double> l, fj, f1j, m, m1, s, s1, s2, fi, f1i, ej, e1j, e2j, n, n1, a, a1, sum1, sum2, cj, c1;
      std::complex<double> k0, k1, k2, c_0, c_1, c_2, c2, n2, d2, d, d1;
      double f, arg, yp, ym, u1, aa1, zcomp1, aa2, u2, zcomp2, zcomp, fac, dfac, pi, pyf, s3;
      vector<double> kr(3), ki(3);
      int i, j;
      pi = 3.14159265359796;
      std::complex<double> c(cos(2*pi/3),sin(2*pi/3));
      std::complex<double> sum(0,0);

      f = 3+3*y-y*y;
      arg=pow((1+2*pow(y,4)/(f*f)),0.5);
      yp= pow((arg+1),(1/3));
      ym=-pow((arg-1),(1/3));

        //for(j=0;j<3;j++)
        //{
         //cj = pow(c,j);
         pyf = pow((2*y*f),(1/3));
         c_0 = pow(c,0);
         c_1 = pow(c,1);
         c_2 = pow(c,2);
         k0 = -2*y+pyf*(yp*c_0+ym/c_0)/(1-y);
         k1 = -2*y+pyf*(yp*c_1+ym/c_1)/(1-y);
         k2 = -2*y+pyf*(yp*c_2+ym/c_2)/(1-y);

         kr[0] = k0.real();
         kr[1] = k1.real();
         kr[2] = k2.real();
         ki[0] = k0.imag();
         ki[1] = k1.imag();
         ki[2] = k2.imag();


         vector< complex<double> > t, tm;
         for(i=0; i<3; i++)
         {
         t.push_back( complex<double>( kr[i], ki[i] ) );
         tm.push_back( complex<double>( -kr[i], -ki[i] ) );
         }

         //std::complex<double>   t[2] = {-2*y+pow((2*y*f),(1/3))*(yp*cj+ym/cj)/(1-y)};
        //}
//=================first order contribution===========================
      for(i=0;i<3;i++)
        {
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         l = (1.0+0.5*y)*t[i]+1.0+2.0*y;
         sum = sum + l*(xlam*t[i]-1.0)*exp((xlam-1.0)*t[i])/(t[i]*s1);
        }

      u1= -1 - 12*y*sum.real();
      aa1=u1;

      sum = sum - sum;
      for(i=0;i<3;i++)
        {
         s1 = 3.0*pow((1-y),2.0)*t[i]*t[i]+12.0*y*(1.0-y)*t[i]+18.0*y*y;
         s2 = 6.0*pow((1.0-y),2.0)*t[i]+12.0*y*(1.0-y);
         a = xlam*(1.0-y)*(1.0+2.0*y)*t[i]*t[i]*t[i]*t[i] +(xlam*pow((1.0+2.0*y),2.0)-(1.0-y)*(1.0+2.0*y))*t[i]*t[i]*t[i] -pow((1.0+2.0*y),2.0)*t[i]*t[i];
         a1 = 4.0*xlam*(1.0-y)*(1.0+2.0*y)*t[i]*t[i]*t[i] +3.0*(xlam*pow((1.0+2.0*y),2.0)-(1.0-y)*(1.0+2.0*y))*t[i]*t[i] -2.0*pow((1.0+2.0*y),2.0)*t[i];
         sum = sum + ((a1+(xlam-1.0)*a)/(s1*s1)-s2*a/(s1*s1*s1))*exp(t[i]*(xlam-1.0));
        }

      zcomp1= -12*y*sum.real();

//==================second order contribution=================

      sum1 = sum1 - sum1;
      sum2 = sum2 - sum2;

      for(i=0;i<3;i++)
        {
         s = pow((1-y),2)*t[i]*t[i]*t[i]+6*y*(1-y)*t[i]*t[i] + 18*y*y*t[i]-12*y*(1+2*y);
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         s2 = 6*pow((1-y),2)*t[i]+12*y*(1-y);
         m = xlam*pow(t[i],5)-pow(t[i],4);
         c1 = 5.0*xlam*pow(t[i],4.0)-4.0*pow(t[i],3.0);
         m1 = c1+(xlam-1)*m;
         fi = m/(s1*s1);
         f1i = m1/(s1*s1)-m*s2/(s1*s1*s1);

         sum1 = sum1+(m*f1i/(s*s)-(m1*s-2.0*m*s1)*fi/(s*s*s));

         for(j=0;j<3;j++)
         {
            s1 = 3*pow((1-y),2)*t[j]*t[j]+12*y*(1-y)*t[j]+18*y*y;
            s2 = 6*pow((1-y),2)*t[j]+12*y*(1-y);
            m = xlam*pow(t[j],5)-pow(t[j],4);
            c1 = 50*xlam*pow(t[j],4.0)-4.0*pow(t[j],3.0);
            m1 = c1+(xlam-1)*m;
            fj = m/(s1*s1);
            f1j = m1/(s1*s1)-m*s2/(s1*s1*s1);

            sum2 = sum2+(f1i*f1j-(fi*f1j+f1i*fj)/(t[i]+t[j]) +2.0*fi*fj/pow((t[i]+t[j]),2.0)) *exp((t[i]+t[j])*(xlam-1))/(t[i]+t[j]);
         }
        }

      u2 = -12*y*pow((1-y),8)*(sum1+sum2).real();
      aa2 = 0.5*u2;
//====================================================
      sum1 = sum1 - sum1;
      sum2 = sum2 - sum2;

      for(i=0;i<3;i++)
        {
         s1 = 3*pow((1-y),2)*t[i]*t[i]+12*y*(1-y)*t[i]+18*y*y;
         s2 = 6*pow((1-y),2)*t[i]+12*y*(1-y);
         m = xlam*pow(t[i],5)-pow(t[i],4);
         c1 = 5.0*xlam*pow(t[i],4.0)-4.0*pow(t[i],3.0);
         m1 = c1+(xlam-1)*m;
         fi = m/(s1*s1);
         f1i = m1/(s1*s1)-m*s2/(s1*s1*s1);

         s = pow((1-y),2)*tm[i]*tm[i]*tm[i]+6*y*(1-y)*tm[i]*tm[i] + 18*y*y*tm[i]-12*y*(1+2*y);
         s1 = 3*pow((1-y),2)*tm[i]*tm[i]+12*y*(1-y)*tm[i]+18*y*y;
         c = xlam*pow(tm[i],5)-pow(tm[i],4);
         d = pow((1.0-y),3.0)*tm[i]*tm[i]*tm[i] -6.0*y*(1.0-y)*(3.0+y)*tm[i]*tm[i] -18.0*y*y*(7.0+y)*tm[i] +12.0*y*(3.0+19.0*y+2.0*y*y);
         n = c*d;
         c1 = 5.0*xlam*pow(tm[i],4.0)-4.0*pow(tm[i],3.0);
         d1 = 3.0*pow((1-y),3.0)*tm[i]*tm[i] -12.0*y*(1.0-y)*(3.0+y)*tm[i] -18.0*y*y*(7.0+y);
         n1 = c1*d+c*d1+(xlam-1)*c*d;

         sum1 = sum1+(n*f1i/(s*s*s)-(n1*s-3.0*n*s1)*fi/(s*s*s*s));

         for(j=0;j<3;j++)
            {
            c = xlam*pow(t[j],5.0)-pow(t[j],4.0);
            d = pow((1-y),3)*t[j]*t[j]*t[j] -6*y*(1-y)*(3+y)*t[j]*t[j] -18*y*y*(7+y)*t[j] +12*y*(3+19*y+2*y*y);
            n = c*d;
            s1 = 3*pow((1-y),2)*t[j]*t[j]+12*y*(1-y)*t[j]+18*y*y;
            ej = n/(s1*s1*s1);

            s2 = 6*pow((1-y),2)*t[j]+12*y*(1-y);
            c1 = 5.0*xlam*pow(t[j],4.0)-4.0*pow(t[j],3.0);
            d1 = 3*pow((1-y),3)*t[j]*t[j] -12*y*(1-y)*(3+y)*t[j] -18*y*y*(7+y);
            n1 = c1*d+c*d1+(xlam-1)*c*d;
            e1j = n1/(s1*s1*s1)-3.0*n*s2/(2.0*pow(s1,4.0));

            c2=20.0*xlam*pow(t[j],3.0)-12.0*pow(t[j],2.0);
            d2=6.0*pow((1-y),3)*t[j] -12*y*(1-y)*(3+y);
            n2 = c2*d+c*d2+pow((xlam-1),2.0)*c*d+2.0*c1*d1 +2.0*(xlam-1.0)*c1*d+2.0*(xlam-1.0)*c*d1;
            s3 = 6.0*pow((1-y),2);
            e2j = n2/(s1*s1*s1)-3.0*n1*s2/pow(s1,4.0) +n*(3.0*s2*s2-s1*s3)/pow(s1,5.0);

            sum2=sum2+(0.5*f1i*e2j -(fi*e2j+2.0*f1i*e1j)/(2.0*(t[i]+t[j])) +(2.0*fi*e1j+f1i*ej)/pow((t[i]+t[j]),2.0) -(3.0*fi*ej)/pow((t[i]+t[j]),3.0)) *exp((t[i]+t[j])*(xlam-1))/(t[i]+t[j]);
            }
        }

      zcomp2 = -6*y*pow((1-y),7)*(sum1+sum2).real();

      zcomp = zcomp1 + zcomp2/tstar;

      fac=(aa1+aa2/tstar)/y;

      dfac=zcomp/(y*y)-fac/y;

      return dfac;
}

double pressref(double sig1, double sig2, double sig12, double rho1, double rho2)
{
      double cnst, cnst1, pi, xi0, xi1, xi2, xi3, d12, g12, sum, term1, term2, term3, bp, termna;
      pi = 3.14159265359796;

      cnst=pi/6;
      cnst1=6/pi;
      xi0=cnst*(rho1+rho2);
      xi1=cnst*(rho1*sig1+rho2*sig2);
      xi2=cnst*(rho1*pow(sig1,2)+rho2*pow(sig2,2));
      xi3=cnst*(rho1*pow(sig1,3)+rho2*pow(sig2,3));

      d12=0.5*(sig1+sig2);
      g12=1/(1-xi3)+1.5*xi2/pow((1-xi3),2)*(sig1*sig2/d12) +0.5*xi2*xi2/pow((1-xi3),3)*pow((sig1*sig2/d12),2);

      sum = xi3/pow((1-xi3),2)+1.5*(sig1*sig2/d12)*xi2*(1+xi3)/pow((1-xi3),3) +0.5*pow((sig1*sig2/d12),2) *xi2*xi2*(2+xi3)/pow((1-xi3),4);

      term1 = xi0/(1-xi3);
      term2 = 3*xi1*xi2/pow(1-xi3,2);
      term3 = (3*pow(xi2,3)-xi3*pow(xi2,3))/pow((1-xi3),3);

      termna = 4*pi*rho1*rho2*d12*g12*(sig12-d12) + 4*pi*rho1*rho2*d12*(sig12-d12)*sum;

      bp = cnst1*(term1+term2+term3) + termna;

      return bp;
}

VectorXd data_msa(int comp)
{
    VectorXd data(7);

    if(comp==1) // CO2
    {
    data[0] = 0.16179E-02; //a111
    data[1] = 0.22172E+01; //a112
    data[2] = 0.20000E+01; //a113
    data[3] = 0.13311E-04; //b111
    data[4] = 0.00000E+00; //b112
    data[5] = 0.00000E+00; //b113
    data[6] = 0.16500E+01; //xlam
    }

    if(comp==2) // n-butane
    {
    data[0] = 0.22613E-02; //a111
    data[1] = 0.36637E+01; //a112
    data[2] = 0.20000E+01; //a113
    data[3] = 0.36177E-04; //b111
    data[4] = 0.00000E+00; //b112
    data[5] = 0.00000E+00; //b113
    data[6] = 0.16500E+01; //xlam
    }

return data;

}

double MSA_pressure(double T, double R, double rho1, double rho2)
{
      double beta, sig1, sig2, sig12, pi, y1, y2, y, t11, t12, t22, cnst, cnst1, xkmix, xlmix;
      double term1, term2, term3, p0, p1, p, bp0, f11, f12, f22, df11, df22, df12;
      double b11, b22, b12, a11, a22, a12, xlam11, xlam12, xlam22, term11, term22, term12;
      double b111, b112, b113, a111, a112, a113, b221, b222, b223, a221, a222, a223;
      VectorXd data(7);

      pi = 3.14159265359796;
      cnst = pi/6;
      cnst1 = 6/pi;

      data = data_msa(1);
      a111 = data[0];
      a112 = data[1];
      a113 = data[2];
      b111 = data[3];
      b112 = data[4];
      b113 = data[5];
      xlam11 = data[6];

      data = data_msa(2);
      a221 = data[0];
      a222 = data[1];
      a223 = data[2];
      b221 = data[3];
      b222 = data[4];
      b223 = data[5];
      xlam22 = data[6];

      xlmix = 0.00000E-03;
      xkmix = 0.18504E+00;

      xlam12 = (xlam11+xlam22)/2;
      b11 = b111 + b112/T + b113/pow(T,2);
      b22 = b221 + b222/T + b223/pow(T,2);
      a11 = (a111+pow((a112/T),a113));
      a22 = (a221+pow((a222/T),a223));

      sig1 = pow((b11/cnst),(1.0/3.0));
      sig2 = pow((b22/cnst),(1.0/3.0));
      sig12 = 0.5*(sig1+sig2)*(1.0-xlmix);
      b12 = cnst*pow(sig12,3.0);
      a12 = pow((a11*a22),0.5)*(1.0-xkmix);

      beta = 1/(R*T);

      sig1 = pow(6*b11/pi,1/3);
      sig2 = pow(6*b22/pi,1/3);
      sig12= pow(6*b12/pi,1/3);

      y1 = b11*rho1;
      y2 = b22*rho2;
      y = y1+y2;

      t11=1/(a11*beta);
      t12=1/(a12*beta);
      t22=1/(a22*beta);

      f11 = script_f(t11,y,xlam11);
      f12 = script_f(t12,y,xlam12);
      f22 = script_f(t22,y,xlam22);

      df11 = script_df(t11,y,xlam11);
      df12 = script_df(t12,y,xlam12);
      df22 = script_df(t22,y,xlam22);

      cout << "rho1 = " << rho1 << endl;
      cout << "rho2 = " << rho2 << endl;
      cout << "y1 = " << y1 << endl;
      cout << "y2 = " << y2 << endl;
      cout << "y = " << y << endl;
      cout << "f11 = " << f11 << endl;
      cout << "f22 = " << f22 << endl;
      cout << "f12 = " << f12 << endl;
      cout << "df11 = " << df11 << endl;
      cout << "df22 = " << df22 << endl;
      cout << "df12 = " << df12 << endl;

      bp0 = pressref(rho1,rho2,sig1,sig2,sig12);

      p0=bp0/beta;

      term11=a11*b11*rho1*rho1*(f11+y*df11);
      term22=a22*b22*rho2*rho2*(f22+y*df22);
      term12=a12*b12*rho1*rho2*(f12+y*df12);
      p1=term11+term22+2*term12;

      p=p0+p1;

      cout << "p0 = " << p0 << endl;
      cout << "p1 = " << p1 << endl;
      cout << "b11 = " << b11 << endl;
      cout << "b22 = " << b22 << endl;
      cout << "a11 = " << a11 << endl;
      cout << "a22 = " << a22 << endl;
      cout << "beta = " << beta << endl;

      return p;
}

MSA_R(double T, double sigma, double eps, double kB)
{
    double c1, c2, c3, R, Tr;

    Tr = kB*T/eps*6.023e23;

    c1 = 1.1287;
    c2 = -0.05536;
    c3 = 0.0007278;

    R = pow(2,1/6) * pow((1+ pow((1+ (Tr+c2*pow(Tr,2)+c3*pow(Tr,4))/c1),1/2)),-1/6) * sigma;
    return R;
}

MSA_S(double t, double eta)
{
    double S;
    S = pow(1-eta,2)*pow(t,3) + 6*eta*(1-eta)*pow(t,2) + 18*pow(eta,2)*t - 12*eta*(1+2*eta);
    return S;
}

MSA_L(double t, double eta)
{
    double L;
    L = (1+eta/2)*t + 1 + 2*eta;
    return L;
}

MSA_Q(double t, double eta)
{
    double Q, L, S;

    S = MSA_S(t,eta);
    L = MSA_L(t,eta);

    Q = (S+12*eta*L*exp(-t))/pow(1-eta,2)/pow(t,3);
    return Q;
}


helmholtz_MSA(double rho, double T, double sigma, double eps, double kB)
{
    double pi, eta, R, a0, a1, a2, a, beta;
    double a11, a12, a13, g0R, Lz1, Lz2, Qz1, Qz2, z1, z2, a21, a22;
    double k0, k1, k2;

    pi = 3.14159265359;
    k0 = 2.1714*sigma;
    z1 = 2.9637/sigma;
    z2 = 14.0167/sigma;

    R = MSA_R(T, sigma, eps, kB);

    eta = pi*rho*pow(R,3)/6;
    beta = 1/kB/T;

    k1 = k0 * exp(z1*(sigma-R));
    k2 = k0 * exp(z2*(sigma-R));

    g0R = (1+eta/2)/pow(1-eta,2);

    Lz1 = MSA_L(z1*R,eta);
    Lz2 = MSA_L(z2*R,eta);
    Qz1 = MSA_Q(z1*R,eta);
    Qz2 = MSA_Q(z2*R,eta);

    a0 = (4*eta-3*pow(eta,2))/pow(1-eta,2);

    a11 = -12*eta*beta*eps/pow(R,3) * (k1*(Lz1/pow(z1,2)/pow(1-eta,2)/Qz1 - (1+z1*R)/pow(z1,2)) -k2*(Lz2/pow(z2,2)/pow(1-eta,2)/Qz2 - (1+z2*R)/pow(z2,2)));
    a12 = +48*eta*beta*eps * (1/9*pow(sigma/R,6) - 1/3*pow(sigma/R,6));
    a13 = -48*eta*beta*eps*g0R * (1/9*pow(sigma/R,12) - 1/3*pow(sigma/R,6) + 2/9*pow(sigma/R,3));
    a1 = a11 + a12 + a13;

    a21 = -6*eta*pow(beta,2)*pow(eps,2)/pow(R,3) * (pow(k1,2)/2/z1/pow(Qz1,4) + pow(k2,2)/2/z2/pow(Qz2,4) - 2*k1*k2/(z1+z2)/pow(Qz1,2)/pow(Qz2,2));
    a22 = -24*eta*pow(beta,2)*pow(eps,2) * (k1/R/pow(Qz1,2) - k2/R/pow(Qz2,2)) * (1/9*pow(sigma/R,12) - 1/3*pow(sigma/R,6) + 2/9*pow(sigma/R,3));
    a2 = a21 + a22;

    a = a0 + a1 + a2;

    return a;
}


#endif // MSA_H_INCLUDED
