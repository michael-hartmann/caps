#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

/**
 * define fresnel_t as either TE or TM.
 */
typedef enum { TE, TM } polarization_t;

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "hash-table.h"
#include "matrix.h"
#include "constants.h"


/* default values */
#define CASIMIR_MINIMUM_LDIM 20 /**< minimum value for lmax */
#define CASIMIR_FACTOR_LDIM 5   /**< by default: lmax=ceil(5/LbyR) */
#define CASIMIR_EPSREL 1e-8  /**< default relative error for integration */

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
    double LbyR; /**< \f$L/R\f$, where \f$R\f$ is the radius of the sphere and \f$L\f$ is the separation of plane and sphere. */
    double y;    /**< log(R/(R+L)/2) */
    /*@}*/

    /**
     * @name reflection properties of the mirrors
     */
     /*@{*/
    double (*epsilonm1)(double xi, void *userdata);
    void *userdata;

    void (*rp)(struct casimir *self, double nT, double k, double *r_TE, double *r_TM);
    void (*lnab)(struct casimir *self, double nT, int l, double *lna, double *lnb);
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


typedef struct {
    casimir_t *casimir;
    int m;
    double tau,epsrel;
    HashTable *hash_table_I;
    HashTable *hash_table_K;
    bool is_pc;
} integration_t;

typedef struct {
    double LbyR, alpha, omegap, epsrel;
    HashTable *cache;
    HashTable *cache_ratio;
} integration_plasma_t;


typedef struct
{
    casimir_t *casimir;
    int m, lmin;
    integration_t *integration;
    integration_plasma_t *integration_plasma;
    double nT;
    double *al, *bl;
} casimir_M_t;


/* prototypes */
int  casimir_compile_info(char *str, size_t size);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

double casimir_lnLambda(int l1, int l2, int m);

casimir_t *casimir_init(double LbyR);
void casimir_free(casimir_t *self);

void casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata);

void casimir_set_rp(casimir_t *self, void (*rp)(struct casimir *self, double nT, double k, double *r_TE, double *r_TM));
void casimir_set_lnab(casimir_t *self, void (*lnab)(struct casimir *self, double nT, int l, double *lna, double *lnb));

int casimir_get_ldim(casimir_t *self);
int casimir_set_ldim(casimir_t *self, int ldim);

detalg_t casimir_get_detalg(casimir_t *self);
void casimir_set_detalg(casimir_t *self, detalg_t detalg);

double casimir_get_epsrel(casimir_t *self);
int    casimir_set_epsrel(casimir_t *self, double epsrel);

void casimir_lnab(casimir_t *self, double nT, int l, double *lna, double *lnb);
void casimir_lnab_perf(casimir_t *self, double nT, int l, double *lna, double *lnb);

double casimir_kernel_M(int i, int j, void *args_);

casimir_M_t *casimir_M_init(casimir_t *self, int m, double nT);
double casimir_M_elem(casimir_M_t *self, int l1, int l2, char p1, char p2);
void casimir_M_free(casimir_M_t *self);
matrix_t *casimir_M(casimir_t *self, double nT, int m);

double casimir_logdetD(casimir_t *self, double nT, int m);

void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM);

int casimir_estimate_lminmax(casimir_t *self, int m, size_t *lmin_p, size_t *lmax_p);

double casimir_epsilonm1(casimir_t *self, double xi);
double casimir_epsilonm1_perf(double xi, void *userdata);
double casimir_epsilonm1_drude(double xi, void *userdata);

double casimir_logdetD0_drude(casimir_t *casimir);
double casimir_logdetD0_pc(casimir_t *casimir, double eps);
double casimir_logdetD0_plasma(casimir_t *casimir, double omegap, double eps);

double casimir_kernel_M0_EE(int i, int j, void *args);
double casimir_kernel_M0_MM(int i, int j, void *args);
double casimir_kernel_M0_MM_plasma(int i, int j, void *args_);

void casimir_logdetD0(casimir_t *self, int m, double omegap, double *EE, double *MM, double *MM_plasma);

#endif
