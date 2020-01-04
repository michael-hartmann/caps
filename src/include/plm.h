#ifndef PLM_H
#define PLM_H

#ifdef __cplusplus
extern "C" {
#endif

#include "attributes.h"

double lnPlm(int l, int m, double x) __attribute__ ((pure));
double lnPlm_upwards(int l, int m, double x) __attribute__ ((pure));
double lnPlm_downwards(int l, int m, double x) __attribute__ ((pure));

double Plm_continued_fraction(const long l, const long m, const double x) __attribute__ ((pure));

double lnPl(int l, double x) __attribute__ ((pure));

double dlnPlm(int l, int m, double x, double *d2lnPlm);

#ifdef __cplusplus
}
#endif

#endif
