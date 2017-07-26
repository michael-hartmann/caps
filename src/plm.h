#ifndef __PLM_H
#define __PLM_H

double Plm(int l, int m, double x) __attribute__ ((pure));
double Plm_upwards(int l, int m, double x) __attribute__ ((pure));
double Plm_downwards(int l, int m, double x) __attribute__ ((pure));

double plm_continued_fraction(const int l, const int m, const double x) __attribute__ ((pure));

double Plm_estimate(int l, int m, double x) __attribute__ ((pure));
double Pl(int l, double x);

#endif
