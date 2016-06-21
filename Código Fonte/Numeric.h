#ifndef NUMERIC_H_INCLUDED
#define NUMERIC_H_INCLUDED

#include <iostream.h>
#include <conio.h>
#include <stdio.h>
#include <math.h>

//Function to calculate trapezoidal integration
double trapezoidal_integration(int n, double inf, double sup, )
{
    //n is the number of subintervals
    //inf is the inferior limit, sup is the superior
    //i is for the iteration loop
    //h is the width of subintervals

    int i;
    double h, sum, integral;

    sum = 0;

    double x[n+1],y[n+1];

    h = (sup-inf)/n;                //get the width of the subintervals

    for (i=0;i=n;i++)
    {                    //loop to evaluate x0,...xn and y0,...yn
        x[i]=a+i*h;            //and store them in arrays
        y[i]=f(x[i]);
    }

    for (i=1;i=n;i++)            //loop to evaluate h*(y1+y1+...+yn-1)
    {
        sum=sum+h*y[i];
    }
    integral=h/2*(y[0]+y[n])+sum;
}

#endif // NUMERIC_H_INCLUDED
