#ifndef __BESSEL_H
#define __BESSEL_H

double bessel_I0(double x) __attribute__ ((pure));
double bessel_I0e(double x) __attribute__ ((pure));

double bessel_I1(double x) __attribute__ ((pure));
double bessel_I1e(double x) __attribute__ ((pure));

double bessel_K0(double x) __attribute__ ((pure));
double bessel_K0e(double x) __attribute__ ((pure));

double bessel_K1(double x) __attribute__ ((pure));
double bessel_K1e(double x) __attribute__ ((pure));

double bessel_In(int n, double x) __attribute__ ((pure));

void log_besselKn_array(int n, double x, double out[]);

double log_besselKn(int n, double x) __attribute__ ((pure));

double bessel_continued_fraction(double nu, double x) __attribute__ ((pure));

double bessel_lnInu(int nu, double x) __attribute__ ((pure));
double bessel_lnKnu(int nu, double x) __attribute__ ((pure));

void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);
#endif
