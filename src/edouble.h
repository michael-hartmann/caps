#ifndef __EDOUBLE_H
#define __EDOUBLE_H

#ifdef EXTENDED_DOUBLE
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
#else
    #if defined(__ICC) || defined(__INTEL_COMPILER)
        #define CASIMIR_ARITHMETICS "icc _Quad"
        #define edouble _Quad

        /* define prototypes. without these prototypes icc will return nan. */
        _Quad __logq(_Quad);
        #define loge __logq

        _Quad __cosq(_Quad);
        #define cose __cosq

        _Quad __sinq(_Quad);
        #define sine __sinq

        _Quad __expq(_Quad);
        #define expe __expq

        _Quad __gammaq(_Quad);
        #define gammae __gammaq

        _Quad __lgammaq(_Quad);
        #define lgammae __lgammaq

        _Quad __sqrtq(_Quad);
        #define sqrte __sqrtq

        _Quad __log1pq(_Quad);
        #define log1pe __log1pq

        _Quad __fabsq(_Quad);
        #define fabse __fabsq

        _Quad __copysignq(_Quad, _Quad);
        #define copysigne __copysignq

    #elif defined(__GNUC__) || defined(__GNUG__)
        #define CASIMIR_ARITHMETICS "gcc __float128"

        #include <quadmath.h>
        #define edouble __float128

        #define loge logq
        #define cose cosq
        #define sine sinq
        #define expe expq
        #define gammae gammaq
        #define lgammae lgammaq
        #define sqrte sqrtq
        #define log1pe log1pq
        #define fabse fabsq
        #define copysigne copysignq
    #else
        #error "I'm sorry, but quad precision is only supported with gcc or icc at the moment."
    #endif
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
    #define COMPILER "icc"
#elif defined(__GNUC__) || defined(__GNUG__)
    #define COMPILER "gcc"
#else
    #define COMPILER "unknown"
#endif

#endif
