#ifndef __FLOATTYPES_H
#define __FLOATTYPES_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double      float64;
typedef long double float80;
typedef __float128  float128;

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

#define log80      logl
#define exp80      expl
#define sqrt80     sqrtl
#define log1p80    log1pl
#define fabs80     fabsl
#define sin80      sinl
#define cos80      cosl
#define gamma80    gammal
#define lgamma80   lgammal
#define copysign80 copysignl

#define log128      logq
#define exp128      expq
#define sqrt128     sqrtq
#define log1p128    log1pq
#define fabs128     fabsq
#define sin128      sinq
#define cos128      cosq
#define gamma128    gammaq
#define lgamma128   lgammaq
#define copysign128 copysignq

#endif
