#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

/**
 * define sign_t as a signed char, because "char can be either signed or
 * unsigned depending on the implementation"
 */
typedef signed char sign_t;

#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>

#include "matrix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define M_LOG2 0.6931471805599453
#define M_LOGPI 1.1447298858494002

/* default values */
#define CASIMIR_PRECISION 1e-12      /**< default precision */
#define CASIMIR_IDLE 1000            /**< idle time in µs */
#define CASIMIR_MINIMUM_LMAX 20      /**< minimum value for lmax */
#define CASIMIR_FACTOR_LMAX 5        /**< by default: lmax=ceil(5/LbyR) */

#ifndef CASIMIR_DETALG
#define CASIMIR_DETALG "LU" /**< default algorithm for matrix decomposition */
#endif


/**
 * Cache for Mie coefficients.
 */
typedef struct
{
    double *ln_al;  /**< list of Mie coefficients \f$a_\ell\f$ (logarithms) */
    sign_t *sign_al; /**< list of signs of Mie coefficients \f$a_\ell\f$ */
    double *ln_bl;  /**< list of Mie coefficients \f$b_\ell\f$ (logarithms) */
    sign_t *sign_bl; /**< list of signs of Mie coefficients \f$b_\ell\f$ */
} casimir_mie_cache_entry_t;

typedef struct {
    int nmax;                            /**< maximum value of n */
    int lmax;                            /**< maximum value of \f$\ell\f$ */
    casimir_mie_cache_entry_t **entries; /**< entries */
    pthread_mutex_t mutex;               /**< mutex for Mie cache */
} casimir_mie_cache_t;

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
    int lmax;               /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    int cores;              /**< number of thread that should be used */
    double precision;       /**< precision \f$\epsilon_p\f$ */
    double threshold;       /**< XXX TBD */
    pthread_t **threads;    /**< list of pthread objects */

    char detalg[128];       /**< algorithm to calculate determinant */

    casimir_mie_cache_t *mie_cache; /**< Mie chache */

    bool verbose; /**< verbose flag */

    /* parameters that you usually do not want to change */
    bool debug; /**< debug flag for more information */

    pthread_mutex_t mutex; /**< mutex for printing */
    /*@}*/
} casimir_t;


/**
 * thread object.
 */
typedef struct
{
    casimir_t *self; /**< pointer to Casimir object */
    int n;           /**< Matsubara term */
    int nmax;        /**< maximum number of n */
    double value;    /**< free energy for Matsubara term n*/
} casimir_thread_t;


/* prototypes */
int  casimir_compile_info(char *str, size_t size);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

int casimir_vfprintf(casimir_t *self, FILE *stream, const char *format, va_list args);
int casimir_debug(casimir_t *self, const char *format, ...);
int casimir_verbose(casimir_t *self, const char *format, ...);

double casimir_lnLambda(int l1, int l2, int m);

int casimir_init(casimir_t *self, double LbyR, double T);
void casimir_free(casimir_t *self);

void casimir_set_debug(casimir_t *self, bool debug);
bool casimir_get_debug(casimir_t *self);

void casimir_set_verbose(casimir_t *self, bool verbose);
bool casimir_get_verbose(casimir_t *self);

int casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata);

int casimir_get_lmax(casimir_t *self);
int casimir_set_lmax(casimir_t *self, int lmax);

int casimir_get_detalg(casimir_t *self, char detalg[128]);
int casimir_set_detalg(casimir_t *self, const char *detalg);

int casimir_get_cores(casimir_t *self);
int casimir_set_cores(casimir_t *self, int cores);

double casimir_get_precision(casimir_t *self);
int    casimir_set_precision(casimir_t *self, double precision);

void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0);
void casimir_lnab(casimir_t *self, int n, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b);
void casimir_lnab_perf(casimir_t *self, int n, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b);

double casimir_F_n(casimir_t *self, const int n, int *mmax);
double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_t *self);
void casimir_mie_cache_alloc(casimir_t *self, int n);
void casimir_mie_cache_clean(casimir_t *self);
void casimir_mie_cache_get(casimir_t *self, int l, int n, double *ln_a, sign_t *sign_a, double *ln_b, sign_t *sign_b);
void casimir_mie_cache_free(casimir_t *self);

void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
matrix_t *casimir_M(casimir_t *self, int n, int m);
double casimir_logdetD(casimir_t *self, int n, int m);

void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM);

double casimir_epsilonm1_perf(double xi, void *userdata);
double casimir_epsilonm1_drude(double xi, void *userdata);

#endif
