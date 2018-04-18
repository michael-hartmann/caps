#ifndef __PLM_H
#define __PLM_H

double lnPlm(int l, int m, double x) __attribute__ ((pure));
double lnPlm_upwards(int l, int m, double x) __attribute__ ((pure));
double lnPlm_downwards(int l, int m, double x) __attribute__ ((pure));

double Plm_continued_fraction(const long l, const long m, const double x) __attribute__ ((pure));

double lnPl(int l, double x);

double dlnPlm(int l, int m, double x, double *d2lnPlm);

#endif
