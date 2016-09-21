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

#include "floattypes.h"

/* default values */
#define CASIMIR_PRECISION 1e-12      /**< default precision */
#define CASIMIR_IDLE 1000            /**< idle time in µs */
#define CASIMIR_MINIMUM_LMAX 20      /**< minimum value for lmax */
#define CASIMIR_FACTOR_LMAX 5        /**< by default: lmax=ceil(5/LbyR) */
#define CASIMIR_TRACE_THRESHOLD 1e-8 /**< threshold to use trace */

#ifndef CASIMIR_DETALG
#define CASIMIR_DETALG "QR_FLOAT64" /**< default algorithm for matrix decomposition */
#endif


/**
 * Cache for Mie coefficients.
 */
typedef struct
{
    float80 *ln_al;  /**< list of Mie coefficients \f$a_\ell\f$ (logarithms) */
    sign_t *sign_al; /**< list of signs of Mie coefficients \f$a_\ell\f$ */
    float80 *ln_bl;  /**< list of Mie coefficients \f$b_\ell\f$ (logarithms) */
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
    double omegap_sphere; /**< plasma frequency \f$\omega_\mathrm{P}\f$ of sphere */
    double omegap_plane;  /**< plasma frequency \f$\omega_\mathrm{P}\f$ of plane */
    double gamma_sphere;  /**< relaxation frequency \f$\gamma\f$ of sphere */
    double gamma_plane;   /**< relaxation frequency \f$\gamma\f$ of plane */
    /*@}*/

    /**
     * @name accuracy and numerical parameters
     */
     /*@{*/
    int lmax;               /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    int cores;              /**< number of thread that should be used */
    double precision;       /**< precision \f$\epsilon_p\f$ */
    double trace_threshold; /**< threshold when Tr M is used as an approximation for log(det(1-M)) */
    pthread_t **threads;    /**< list of pthread objects */

    char detalg[128];       /**< algorithm to calculate determinant */

    casimir_mie_cache_t *mie_cache; /**< Mie chache */

    bool verbose; /**< verbose flag */

    /* parameters that you usually do not want to change */
    bool check_elems;  /**< check if matrix elements of M are finite */
    bool debug;        /**< debug flag for more information */

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


typedef struct
{
    float80 lnA_TE;  /**< logarithm of integral A_TE */
    float80 lnA_TM;  /**< logarithm of integral A_TM */
    float80 lnB_TE;  /**< logarithm of integral B_TE */
    float80 lnB_TM;  /**< logarithm of integral B_TM */
    float80 lnC_TE;  /**< logarithm of integral C_TE */
    float80 lnC_TM;  /**< logarithm of integral C_TM */
    float80 lnD_TE;  /**< logarithm of integral D_TE */
    float80 lnD_TM;  /**< logarithm of integral D_TM */
    sign_t signA_TE; /**< sign of lnA_TE */
    sign_t signA_TM; /**< sign of lnA_TM */
    sign_t signB_TE; /**< sign of lnB_TE */
    sign_t signB_TM; /**< sign of lnB_TM */
    sign_t signC_TE; /**< sign of lnC_TE */
    sign_t signC_TM; /**< sign of lnC_TM */
    sign_t signD_TE; /**< sign of lnD_TE */
    sign_t signD_TM; /**< sign of lnD_TM */
} casimir_integrals_t;

#include <matrix.h>

/* prototypes */
int  casimir_compile_info(char *str, size_t size);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

int casimir_vfprintf(casimir_t *self, FILE *stream, const char *format, va_list args);
int casimir_debug(casimir_t *self, const char *format, ...);
int casimir_verbose(casimir_t *self, const char *format, ...);

float80 casimir_epsilonm1(float80 xi, double omegap, double gamma_);

float80 casimir_lnLambda(int l1, int l2, int m, sign_t *sign);

int casimir_init(casimir_t *self, double LbyR, double T);
void casimir_free(casimir_t *self);

void casimir_set_debug(casimir_t *self, bool debug);
bool casimir_get_debug(casimir_t *self);

void casimir_set_verbose(casimir_t *self, bool verbose);
bool casimir_get_verbose(casimir_t *self);

int casimir_set_drude(casimir_t *self, double omegap_plane, double gamma_plane, double omegap_sphere, double gamma_sphere);
int casimir_get_drude(casimir_t *self, double *omegap_plane, double *gamma_plane, double *omegap_sphere, double *gamma_sphere);

int casimir_get_lmax(casimir_t *self);
int casimir_set_lmax(casimir_t *self, int lmax);

int casimir_get_detalg(casimir_t *self, char detalg[128]);
int casimir_set_detalg(casimir_t *self, const char *detalg);

int casimir_get_cores(casimir_t *self);
int casimir_set_cores(casimir_t *self, int cores);

double casimir_get_precision(casimir_t *self);
int    casimir_set_precision(casimir_t *self, double precision);


double casimir_get_trace_threshold(casimir_t *self);
int    casimir_set_trace_threshold(casimir_t *self, double threshold);

void casimir_lnab0(int l, float80 *a0, sign_t *sign_a0, float80 *b0, sign_t *sign_b0);
void casimir_lnab(casimir_t *self, int n, int l, float80 *lna, float80 *lnb, sign_t *sign_a, sign_t *sign_b);
void casimir_lnab_perf(casimir_t *self, int n, int l, float80 *lna, float80 *lnb, sign_t *sign_a, sign_t *sign_b);

double casimir_F_n(casimir_t *self, const int n, int *mmax);
double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_t *self);
void casimir_mie_cache_alloc(casimir_t *self, int n);
void casimir_mie_cache_clean(casimir_t *self);
void casimir_mie_cache_get(casimir_t *self, int l, int n, float80 *ln_a, sign_t *sign_a, float80 *ln_b, sign_t *sign_b);
void casimir_mie_cache_free(casimir_t *self);

void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
matrix_t *casimir_M(casimir_t *self, int n, int m);
double casimir_logdetD(casimir_t *self, int n, int m);
double casimir_trM(casimir_t *self, int n, int m, void *int_perf);

void casimir_rp(casimir_t *self, float80 nT, float80 k, float80 *r_TE, float80 *r_TM);

#endif
