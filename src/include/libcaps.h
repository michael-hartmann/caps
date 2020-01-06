/**
 * @file   libcaps.h
 * @author Michael Hartmann <caps@speicherleck.de>
 * @date   February, 2019
 */

#ifndef LIBCAPS_H
#define LIBCAPS_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * typoe for polarization: either TE or TM.
 */
typedef enum { TE, TM } polarization_t;

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "attributes.h"
#include "cache.h"
#include "matrix.h"


/* default values */
#define CAPS_MINIMUM_LDIM 20 /**< minimum value for lmax */
#define CAPS_FACTOR_LDIM 5   /**< by default: lmax=ceil(5/LbyR) */
#define CAPS_EPSREL 1e-8  /**< default relative error for integration */

#define CAPS_CACHE_ELEMS 10000000 /**< default number of elems of the cache for I integrals */

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
#define POW_2(x) ((x)*(x))

#define CAPS_PI      3.14159265358979323846 /**< value for \f$\pi=3.141592...\f$ */
#define CAPS_LOG2    0.6931471805599453 /**< \f$\log(2)\f$ */
#define CAPS_LOGPI   1.1447298858494002 /**< \f$\log(\pi)\f$ */
#define CAPS_HBAR    1.0545718e-34   /**< reduced Planck constant \f$\hbar\f$ [m² kg / s] */
#define CAPS_HBAR_EV 6.582119514e-16 /**< reduced Planck constant \f$\hbar\f$ [eV s/rad] */
#define CAPS_KB      1.38064852e-23  /**< Boltzman constant \f$k_\mathrm{B}\f$ [m² kg / ( K s² )] */
#define CAPS_C       299792458.      /**< speed of light \f$c\f$ in vacuum [m/s] */

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

/**
 * The CaPS object. This structure stores all essential information about
 * temperature, geometry and the reflection properties of the mirrors.
 *
 * Do not modify the attributes of the structure yourself!
 */
typedef struct caps
{
    /**
     * @name geometry
     */
     /*@{*/
    double L;    /**< separation of plane and sphere*/
    double R;    /**< radius of sphere */
    double calL; /**< \f$L+R\f$ */
    double LbyR; /**< \f$L/R\f$ */
    double y;    /**< log(R/(R+L)/2) */
    /*@}*/

    /**
     * @name dielectric function of the plate
     */
     /*@{*/
    double (*epsilonm1_plate)(double xi_, void *userdata);
    void *userdata_plate;
    /*@}*/

    /**
     * @name dielectric function of the sphere
     */
     /*@{*/
    double (*epsilonm1_sphere)(double xi_, void *userdata);
    void *userdata_sphere;
    /*@}*/

    /**
     * @name accuracy and numerical parameters
     */
     /*@{*/
    int ldim;        /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    double epsrel;   /**< relative error for integration */
    detalg_t detalg; /**< algorithm to calculate determinant */
    /*@}*/
} caps_t;


/**
 * object for integration over k in matrix elements of round-trip operator
 */
typedef struct {
    caps_t *caps;
    int m;
    double alpha,epsrel;
    cache_t *cache_I;
    double *cache_K[2];
    size_t elems_cache_K;
    bool is_pr;
} integration_t;

typedef struct {
    double LbyR, alpha, omegap_, epsrel;
    cache_t *cache;
    cache_t *cache_ratio;
} integration_plasma_t;


typedef struct
{
    caps_t *caps;
    int m, lmin, ldim;
    integration_t *integration;
    integration_plasma_t *integration_plasma;
    double xi_;
    double *al, *bl;
} caps_M_t;


/* prototypes */
void caps_build(FILE *stream, const char *prefix);
void caps_info(caps_t *self, FILE *stream, const char *prefix);

double caps_lnLambda(int l1, int l2, int m);

caps_t *caps_init(double R, double L);
void caps_free(caps_t *self);

double caps_epsilonm1_plate(caps_t *self, double xi_);
double caps_epsilonm1_sphere(caps_t *self, double xi_);

void caps_set_epsilonm1(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);
void caps_set_epsilonm1_plate(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);
void caps_set_epsilonm1_sphere(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);

int caps_get_ldim(caps_t *self);
int caps_set_ldim(caps_t *self, int ldim);

int caps_detalg_from_string(const char *str, detalg_t *detalg);
int caps_detalg_to_string(detalg_t detalg, const char **str);
detalg_t caps_get_detalg(caps_t *self);
void caps_set_detalg(caps_t *self, detalg_t detalg);

double caps_get_epsrel(caps_t *self);
int    caps_set_epsrel(caps_t *self, double epsrel);

void caps_mie(caps_t *self, double xi_, int l, double *lna, double *lnb);
void caps_mie_perf(caps_t *self, double xi_, int l, double *lna, double *lnb);

double caps_kernel_M(int i, int j, void *args_);

caps_M_t *caps_M_init(caps_t *self, int m, double xi_);
double caps_M_elem(caps_M_t *self, int l1, int l2, char p1, char p2);
void caps_M_free(caps_M_t *self);

double caps_logdetD(caps_t *self, double xi_, int m);

void caps_fresnel(caps_t *self, double xi_, double k, double *r_TE, double *r_TM);

int caps_estimate_lminmax(caps_t *self, int m, size_t *lmin_p, size_t *lmax_p);

double caps_epsilonm1_perf(__attribute__((unused)) double xi_, __attribute__((unused)) void *userdata);
double caps_epsilonm1_drude(double xi_, void *userdata);

double caps_ht_drude(caps_t *caps);
double caps_ht_perf(caps_t *caps, double eps);
double caps_ht_plasma(caps_t *caps, double omegap_, double eps);

double caps_kernel_M0_EE(int i, int j, void *args);
double caps_kernel_M0_MM(int i, int j, void *args);
double caps_kernel_M0_MM_plasma(int i, int j, void *args_);

void caps_logdetD0(caps_t *self, int m, double omegap, double *EE, double *MM, double *MM_plasma);

#ifdef __cplusplus
}
#endif

#endif
