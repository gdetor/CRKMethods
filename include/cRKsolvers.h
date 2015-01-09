#ifndef CRKSOLVERS_H
#define CRKSOLVERS_H
#include <stdlib.h>

double f(double , double );
double fy(double , double );

void eulerStep( double (*)(double,double), double, double, double, int, 
                   double *, double *);
    
void rk4Step( double (*)(double,double), double, double, double, int, double *,
              double *);

void rk4RefinedStep(double (*)(double,double), double (*)(double,double), 
                      double, double, double, int, double *, double *);

void rkFehlbergStep( double (*)(double,double), double (*)(double,double), 
                        double, double, double, int, double *, double *);

#endif /* CRKSOLVERS_H */
