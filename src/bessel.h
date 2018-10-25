#ifndef __BESSEL_H
#define __BESSEL_H

double bessel_I0(double x) __attribute__ ((pure));
double bessel_I1(double x) __attribute__ ((pure));

double bessel_K0(double x) __attribute__ ((pure));
double bessel_K1(double x) __attribute__ ((pure));

double bessel_logI0(double x) __attribute__ ((pure));
double bessel_logK0(double x) __attribute__ ((pure));

double bessel_logI1(double x) __attribute__ ((pure));
double bessel_logK1(double x) __attribute__ ((pure));

double bessel_In(int n, double x) __attribute__ ((pure));
double bessel_Kn(int n, double x) __attribute__ ((pure));

void bessel_logKn_recursive(int nmax, double x, double *logKn);

double bessel_logIn(int n, double x) __attribute__ ((pure));
double bessel_logKn(int n, double x) __attribute__ ((pure));

double bessel_ratioI(double nu, double x) __attribute__ ((pure));

double bessel_logInu_series(double nu, double x) __attribute__ ((pure));
double bessel_logInu_asymp(double nu, double x) __attribute__ ((pure));
double bessel_logKnu_asymp(double nu, double x) __attribute__ ((pure));

void bessel_logInuKnu_half(int nu, const double x, double *logInu_p, double *logKnu_p);
double bessel_logInu_half(int nu, double x) __attribute__ ((pure));
double bessel_logKnu_half(int nu, double x) __attribute__ ((pure));

#endif
