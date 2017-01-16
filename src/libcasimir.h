#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

/**
 * define fresnel_t as either TE or TM.
 */
typedef enum { TE, TM } polarization_t;

#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "matrix.h"

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

/* default values */
#define CASIMIR_PRECISION 1e-12      /**< default precision */
#define CASIMIR_MINIMUM_LMAX 20      /**< minimum value for lmax */
#define CASIMIR_FACTOR_LMAX 5        /**< by default: lmax=ceil(5/LbyR) */

/**
 * The Casimir object. This structure stores all essential information about
 * temperature, geometry and the reflection properties of the mirrors.
 *
 * Do not modify the attributes of the structure yourself!
 */
typedef struct
{
    /**
     * @name geometry and temperature
     */
     /*@{*/
    double RbyScriptL; /**< \f$R/\mathcal{L}\f$, where \f$R\f$ is the radius of the sphere and \f$L\f$ is the separation of plane and sphere. */
    double LbyR;       /**< \f$L/R\f$, where \f$R\f$ is the radius of the sphere and \f$L\f$ is the separation of plane and sphere. */
    double T;          /**< temperature */
    /*@}*/

    /**
     * @name reflection properties of the mirrors
     */
     /*@{*/
    double (*epsilonm1)(double xi, void *userdata);
    void *userdata;
    /*@}*/

    /**
     * @name accuracy and numerical parameters
     */
     /*@{*/
    int ldim;               /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    double precision;       /**< precision \f$\epsilon_p\f$ */
    double threshold;       /**< XXX TBD */

    detalg_t detalg;        /**< algorithm to calculate determinant */

    bool verbose; /**< verbose flag */

    /* parameters that you usually do not want to change */
    bool debug; /**< debug flag for more information */

    pthread_mutex_t mutex; /**< mutex for printing */
    /*@}*/
} casimir_t;


/* prototypes */
int  casimir_compile_info(char *str, size_t size);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

int casimir_vfprintf(casimir_t *self, FILE *stream, const char *format, va_list args);
int casimir_debug(casimir_t *self, const char *format, ...);
int casimir_verbose(casimir_t *self, const char *format, ...);

double casimir_lnLambda(int l1, int l2, int m);

casimir_t *casimir_init(double LbyR);
void casimir_free(casimir_t *self);

void casimir_set_debug(casimir_t *self, bool debug);
bool casimir_get_debug(casimir_t *self);

void casimir_set_verbose(casimir_t *self, bool verbose);
bool casimir_get_verbose(casimir_t *self);

int casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata);

int casimir_get_ldim(casimir_t *self);
int casimir_set_ldim(casimir_t *self, int ldim);

detalg_t casimir_get_detalg(casimir_t *self);
int casimir_set_detalg(casimir_t *self, detalg_t detalg);

double casimir_get_precision(casimir_t *self);
int    casimir_set_precision(casimir_t *self, double precision);

double casimir_get_threshold(casimir_t *self);
int    casimir_set_threshold(casimir_t *self, double threshold);

void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0);
void casimir_lnab(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b);
void casimir_lnab_perf(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b);

void casimir_M0(casimir_t *self, int m, matrix_t **EE, matrix_t **MM);
void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
matrix_t *casimir_M(casimir_t *self, double nT, int m);
double casimir_logdetD(casimir_t *self, double nT, int m);

void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM);

int casimir_estimate_lminmax(casimir_t *self, int m, size_t *lmin_p, size_t *lmax_p);

double casimir_epsilonm1(casimir_t *self, double xi);
double casimir_epsilonm1_perf(double xi, void *userdata);
double casimir_epsilonm1_drude(double xi, void *userdata);

#endif
