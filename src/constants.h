#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#define pow_2(x) ((x)*(x))
#define pow_3(x) ((x)*(x)*(x))
#define pow_4(x) ((x)*(x)*(x)*(x))

/* calculate pow(-1,a) = -1**a */
#define MPOW(a) (1-2*((signed char)(a) & 1))

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

#define SGN(val) ((0 < (val)) - ((val) < 0))

#define pow_2(x) ((x)*(x))
#define pow_3(x) ((x)*(x)*(x))
#define pow_4(x) ((x)*(x)*(x)*(x))

/* calculate pow(-1,a) = -1**a */
#define MPOW(a) (1-2*((signed char)(a) & 1))

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

#define SGN(val) ((0 < (val)) - ((val) < 0))

#ifndef M_PI
#define M_PI 3.14159265358979323846 /**< π */
#endif

#ifndef M_LOG2
#define M_LOG2 0.6931471805599453 /**< log(2) */
#endif

#ifndef M_LOGPI
#define M_LOGPI 1.1447298858494002 /**< log(π) */
#endif

#ifndef M_GM
#define M_GM 1.618033988749895 /**< golden mean, (1+√5)/2 */
#endif

#define CASIMIR_hbar    1.0545718e-34   /**< reduced Planck constant [m² kg / s] */
#define CASIMIR_hbar_eV 6.582119514e-16 /**< reduced Planck constant [eV s] */
#define CASIMIR_kB      1.38064852e-23  /**< Boltzman constant [m² kg / ( K s² )] */
#define CASIMIR_c       299792458.      /**< speed of light [m/s] */

#endif
