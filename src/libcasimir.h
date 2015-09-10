#ifndef __LIBCASIMIR_H
#define __LIBCASIMIR_H

#include <pthread.h>
#include <stdio.h>

#include "edouble.h"

#define HBARC 3.161526510740123e-26 /**< \f$\hbar \mathrm{c}\f$ */
#define KB    1.3806488e-23         /**< \f$k_\mathrm{B}\f$ */

#define CASIMIR_DEFAULT_PRECISION 1e-12
#define CASIMIR_IDLE 1000
#define CASIMIR_FACTOR_LMAX 5

/** define sign_t as a signed char */
typedef char sign_t;


/**
 * Cache for Mie coefficients.
 */
typedef struct
{
    double *ln_al;   /**< list of Mie coefficients \f$a_\ell\f$ (logarithms) */
    sign_t *sign_al; /**< list of signs of Mie coefficients \f$a_\ell\f$ */
    double *ln_bl;   /**< list of Mie coefficients \f$b_\ell\f$ (logarithms) */
    sign_t *sign_bl; /**< list of signs of Mie coefficients \f$b_\ell\f$ */
} casimir_mie_cache_entry_t;

typedef struct {
    int nmax;  /**< maximum value of n */
    int lmax;  /**< maximum value of \f$\ell\f$ */
    casimir_mie_cache_entry_t **entries; /**< entries */
    pthread_mutex_t mutex; /**< mutex for Mie cache */
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
    int integration;     /**< 0: use perfect reflectors, >0: order of Gauss-Laguerre integration */
    int lmax;            /**< truncation value for vector space \f$\ell_\mathrm{max}\f$ */
    int cores;           /**< number of thread that should be used */
    double precision;    /**< precision \f$\epsilon_p\f$ */
    pthread_t **threads; /**< list of pthread objects */

    casimir_mie_cache_t *mie_cache;
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
    edouble lnA_TE;  /**< logarithm of integral A_TE */
    edouble lnA_TM;  /**< logarithm of integral A_TM */
    edouble lnB_TE;  /**< logarithm of integral B_TE */
    edouble lnB_TM;  /**< logarithm of integral B_TM */
    edouble lnC_TE;  /**< logarithm of integral C_TE */
    edouble lnC_TM;  /**< logarithm of integral C_TM */
    edouble lnD_TE;  /**< logarithm of integral D_TE */
    edouble lnD_TM;  /**< logarithm of integral D_TM */
    sign_t signA_TE; /**< sign of lnA_TE */
    sign_t signA_TM; /**< sign of lnA_TM */
    sign_t signB_TE; /**< sign of lnB_TE */
    sign_t signB_TM; /**< sign of lnB_TM */
    sign_t signC_TE; /**< sign of lnC_TE */
    sign_t signC_TM; /**< sign of lnC_TM */
    sign_t signD_TE; /**< sign of lnD_TE */
    sign_t signD_TM; /**< sign of lnD_TM */
} casimir_integrals_t;


/* prototypes */
const char *casimir_compile_info(void);
void casimir_info(casimir_t *self, FILE *stream, const char *prefix);

double casimir_epsilon(double xi, double omegap, double gamma_);
double casimir_lnepsilon(double xi, double omegap, double gamma_);

edouble casimir_lnLambda(int l1, int l2, int m, sign_t *sign);
edouble casimir_lnXi(int l1, int l2, int m, sign_t *sign);

double casimir_F_SI_to_scaled(double F_SI, double ScriptL_SI);
double casimir_F_scaled_to_SI(double F, double ScriptL_SI);
double casimir_T_SI_to_scaled(double T_SI, double ScriptL_SI);
double casimir_T_scaled_to_SI(double T, double ScriptL_SI);

int casimir_init(casimir_t *self, double RbyScriptL, double T);
void casimir_free(casimir_t *self);

int casimir_set_omegap_sphere(casimir_t *self, double omegap);
int casimir_set_omegap_plane(casimir_t *self, double omegap);

double casimir_get_omegap_sphere(casimir_t *self);
double casimir_get_omegap_plane(casimir_t *self);

int casimir_set_gamma_sphere(casimir_t *self, double gamma_);
int casimir_set_gamma_plane(casimir_t *self, double gamma_);

double casimir_get_gamma_sphere(casimir_t *self);
double casimir_get_gamma_plane(casimir_t *self);

void casimir_set_integration(casimir_t *self, int integration);
int casimir_get_integration(casimir_t *self);

int casimir_get_lmax(casimir_t *self);
int casimir_set_lmax(casimir_t *self, int lmax);

int casimir_get_cores(casimir_t *self);
int casimir_set_cores(casimir_t *self, int cores);

double casimir_get_precision(casimir_t *self);
int    casimir_set_precision(casimir_t *self, double precision);

void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0);
void casimir_lnab(casimir_t *self, const int n, const int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b);
double casimir_lna_perf(casimir_t *self, const int l, const int n, sign_t *sign);
double casimir_lnb_perf(casimir_t *self, const int l, const int n, sign_t *sign);

double casimir_F_n(casimir_t *self, const int n, int *mmax);
double casimir_F(casimir_t *self, int *nmax);

void casimir_mie_cache_init(casimir_t *self);
void casimir_mie_cache_alloc(casimir_t *self, int n);
void casimir_mie_cache_clean(casimir_t *self);
void casimir_mie_cache_get(casimir_t *self, int l, int n, double *ln_a, sign_t *sign_a, double *ln_b, sign_t *sign_b);
void casimir_mie_cache_free(casimir_t *self);

void casimir_logdetD0(casimir_t *self, int m, double *EE, double *MM);
double casimir_logdetD(casimir_t *self, int n, int m, void * obj);
double casimir_trM(casimir_t *self, int n, int m, void *integration_obj);

void casimir_rp(casimir_t *self, edouble nT, edouble k, edouble *r_TE, edouble *r_TM);

#endif
