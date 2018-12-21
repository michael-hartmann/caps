/**
 * @file   libcasimir.h
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   December, 2017
 */

#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

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

#include "cache.h"
#include "matrix.h"
#include "constants.h"


/* default values */
#define CASIMIR_MINIMUM_LDIM 20 /**< minimum value for lmax */
#define CASIMIR_FACTOR_LDIM 5   /**< by default: lmax=ceil(5/LbyR) */
#define CASIMIR_EPSREL 1e-8  /**< default relative error for integration */

#define CASIMIR_CACHE_ELEMS 1000000 /**< default number of elems of the cache for I integrals */

/**
 * The Casimir object. This structure stores all essential information about
 * temperature, geometry and the reflection properties of the mirrors.
 *
 * Do not modify the attributes of the structure yourself!
 */
typedef struct casimir
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
     * @name dielectric function of the plate
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
} casimir_t;


/**
 * object for integration over k in matrix elements of round-trip operator
 */
typedef struct {
    casimir_t *casimir;
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
    casimir_t *casimir;
    int m, lmin;
    integration_t *integration;
    integration_plasma_t *integration_plasma;
    double xi_;
    double *al, *bl;
} casimir_M_t;


/* prototypes */
void casimir_build(FILE *stream, const char *prefix);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

double casimir_lnLambda(int l1, int l2, int m);

casimir_t *casimir_init(double R, double L);
void casimir_free(casimir_t *self);

double casimir_epsilonm1_plate(casimir_t *self, double xi_);
double casimir_epsilonm1_sphere(casimir_t *self, double xi_);

void casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);
void casimir_set_epsilonm1_plate(casimir_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);
void casimir_set_epsilonm1_sphere(casimir_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata);

int casimir_get_ldim(casimir_t *self);
int casimir_set_ldim(casimir_t *self, int ldim);

detalg_t casimir_get_detalg(casimir_t *self);
int casimir_set_detalg(casimir_t *self, detalg_t detalg);

double casimir_get_epsrel(casimir_t *self);
int    casimir_set_epsrel(casimir_t *self, double epsrel);

void casimir_mie(casimir_t *self, double xi_, int l, double *lna, double *lnb);
void casimir_mie_perf(casimir_t *self, double xi_, int l, double *lna, double *lnb);

double casimir_kernel_M(int i, int j, void *args_);

casimir_M_t *casimir_M_init(casimir_t *self, int m, double xi_);
double casimir_M_elem(casimir_M_t *self, int l1, int l2, char p1, char p2);
void casimir_M_free(casimir_M_t *self);

double casimir_logdetD(casimir_t *self, double xi_, int m);

void casimir_fresnel(casimir_t *self, double xi_, double k, double *r_TE, double *r_TM);

int casimir_estimate_lminmax(casimir_t *self, int m, size_t *lmin_p, size_t *lmax_p);

double casimir_epsilonm1_perf(__attribute__((unused)) double xi_, __attribute__((unused)) void *userdata);
double casimir_epsilonm1_drude(double xi_, void *userdata);

double casimir_ht_drude(casimir_t *casimir);
double casimir_ht_perf(casimir_t *casimir, double eps);
double casimir_ht_plasma(casimir_t *casimir, double omegap_, double eps);

double casimir_kernel_M0_EE(int i, int j, void *args);
double casimir_kernel_M0_MM(int i, int j, void *args);
double casimir_kernel_M0_MM_plasma(int i, int j, void *args_);

void casimir_logdetD0(casimir_t *self, int m, double omegap, double *EE, double *MM, double *MM_plasma);

#ifdef __cplusplus
}
#endif

#endif
