#ifndef __EDOUBLE_H
#define __EDOUBLE_H

#define CASIMIR_ARITHMETICS "long double"
#define edouble long double

#define loge      logl
#define expe      expl
#define sqrte     sqrtl
#define log1pe    log1pl
#define fabse     fabsl
#define sine      sinl
#define cose      cosl
#define gammae    gammal
#define lgammae   lgammal
#define copysigne copysignl

#if defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#else
    #define COMPILER "unknown"
#endif

#endif
