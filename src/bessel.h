#ifndef __BESSEL_H
#define __BESSEL_H

double besselI0(double x);
double besselI0e(double x);

double besselI1(double x) __attribute__ ((pure));;
double besselI1e(double x) __attribute__ ((pure));;

double besselI(int n, double x) __attribute__ ((pure));;

double bessel_continued_fraction(int nu, double x) __attribute__ ((pure));;

double bessel_lnInu(int nu, double x) __attribute__ ((pure));;
double bessel_lnKnu(int nu, double x) __attribute__ ((pure));;

void bessel_lnInuKnu(int nu, const double x, double *lnInu_p, double *lnKnu_p);
#endif
