#ifndef FCQS_H
#define FCQS_H

#ifdef __cplusplus
extern "C" {
#endif

double fcqs_semiinf(double f(double, void *), void *args, double *epsrel, int *neval, double L, int *ier);
double fcqs_finite(double f(double, void *), void *args, double a, double b, double *epsrel, int *neval, int *ier);

#ifdef __cplusplus
}
#endif

#endif
