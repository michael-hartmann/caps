#ifndef FCQS_H
#define FCQS_H

double fcgs_semiinf(double f(double, void *), void *args, double *epsrel, int *neval, double L, int *ier);
double fcgs_finite(double f(double, void *), void *args, double a, double b, double *epsrel, int *neval, int *ier);

#endif
