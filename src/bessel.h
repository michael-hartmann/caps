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

void log_besselKn_array(int n, double x, double out[]);

double log_besselKn(int n, double x) __attribute__ ((pure));

double bessel_continued_fraction(double nu, double x) __attribute__ ((pure));

double bessel_logInu(int nu, double x) __attribute__ ((pure));
double bessel_logKnu(int nu, double x) __attribute__ ((pure));

void bessel_logInuKnu(int nu, const double x, double *logInu_p, double *logKnu_p);
#endif
