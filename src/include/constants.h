/**
 * @file   constants.h
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   October, 2018
 * @brief  define macros and constants
 */

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

/*! macro to get minimum of two numbers */
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

/*! macro to get maximum of two numbers */
#ifndef MAX
#define MAX(a,b) ((((a))>((b)))?((a)):((b)))
#endif

/*! macro to get sign of numbers */
#define SGN(val) ((0 < (val)) - ((val) < 0))

/*! compute \f$x^2\f$ */
#define pow_2(x) ((x)*(x))

#ifndef M_PI
#define M_PI 3.14159265358979323846 /**< value for \f$\pi=3.141592...\f$ */
#endif

#ifndef M_LOG2
#define M_LOG2 0.6931471805599453 /**< \f$log(2)\f$ */
#endif

#ifndef M_LOGPI
#define M_LOGPI 1.1447298858494002 /**< \f$log(\pi)\f$ */
#endif

#define CASIMIR_hbar    1.0545718e-34   /**< reduced Planck constant \f$\hbar\f$ [m² kg / s] */
#define CASIMIR_hbar_eV 6.582119514e-16 /**< reduced Planck constant \f$\hbar\f$ [eV s/rad] */
#define CASIMIR_kB      1.38064852e-23  /**< Boltzman constant \f$k_B\f$ [m² kg / ( K s² )] */
#define CASIMIR_c       299792458.      /**< speed of light \f$c\f$ in vacuum [m/s] */

#endif
