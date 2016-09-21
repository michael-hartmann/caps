#ifndef __FLOATTYPES_H
#define __FLOATTYPES_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define M_LOG2 0.6931471805599453
#define M_LOGPI 1.1447298858494002


typedef double     float64;
#define log64      log
#define exp64      exp
#define sqrt64     sqrt
#define log1p64    log1p
#define fabs64     fabs
#define sin64      sin
#define cos64      cos
#define gamma64    gamma
#define lgamma64   lgamma
#define copysign64 copysign
#define fmax64     fmax
#define fmin64     fmin
#define pow64      pow

#endif
