/**
 * @file   libcasimir.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */


#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

#include "edouble.h"
#include "integration_drude.h"
#include "integration_perf.h"
#include "libcasimir.h"
#include "matrix.h"
#include "sfunc.h"
#include "utils.h"

static char CASIMIR_COMPILE_INFO[4096] = { 0 };

/** @brief Return string with information about the binary
 *
 * The returned string contains date and time of compilation, the compiler and
 * kind of arithmetics the binary uses.
 *
 * Do not modify or free this string!
 *
 * This function is not thread-safe.
 *
 * @retval description constant string
 */
const char *casimir_compile_info(void)
{
    #ifdef USE_LAPACK
        const char *lapack = "lapack support";
    #else
        const char *lapack = "no lapack support";
    #endif

    snprintf(CASIMIR_COMPILE_INFO, sizeof(CASIMIR_COMPILE_INFO)/sizeof(char),
             "Compiled on %s at %s with %s, using %s, %s",
              __DATE__, __TIME__, COMPILER, CASIMIR_ARITHMETICS, lapack
            );
    return CASIMIR_COMPILE_INFO;
}


/** @brief Print object information to stream
 *
 * This function will print information about the object self to stream.
 * Information include all parameters like \f$R/L\f$, \f$\omega_\mathrm{P}\f$
 * and \f$\gamma\f$ of sphere and plane, as well as maximum value of
 * \f$\ell\f$, precision, number of cores...
 *
 * This function is thread-safe. However, do not modify parameters (e.g. lmax,
 * dielectric properties of plane and sphere...) while calling this function.
 *
 * @param self Casimir object
 * @param stream where to print the string
 * @param prefix if prefix != NULL: start every line with the string contained
 * in prefix
 */
void casimir_info(casimir_t *self, FILE *stream, const char *prefix)
{
    if(prefix == NULL)
        prefix = "";

    fprintf(stream, "%sR/(R+L) = %g\n", prefix, self->RbyScriptL);
    fprintf(stream, "%sL/R     = %g\n", prefix, self->LbyR);
    fprintf(stream, "%sR/L     = %g\n", prefix, 1/self->LbyR);
    fprintf(stream, "%sT       = %g\n", prefix, self->T);

    fprintf(stream, "%somegap_sphere = %g\n", prefix, self->omegap_sphere);
    fprintf(stream, "%somegap_plane  = %g\n", prefix, self->omegap_plane);
    fprintf(stream, "%sgamma_sphere  = %g\n", prefix, self->gamma_sphere);
    fprintf(stream, "%sgamma_plane   = %g\n", prefix, self->gamma_plane);
    fprintf(stream, "%sintegration   = ", prefix);
    if(self->integration <= 0)
        fprintf(stream, "analytic (perfect mirrors)\n");
    else
        fprintf(stream, "%d\n", self->integration);

    fprintf(stream, "%slmax      = %d\n", prefix, self->lmax);
    fprintf(stream, "%sverbose   = %d\n", prefix, self->verbose);
    fprintf(stream, "%scores     = %d\n", prefix, self->cores);
    fprintf(stream, "%sprecision = %g\n", prefix, self->precision);
}


/**
 * @brief Calculate logarithm and sign of prefactor \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of the prefactor for given
 * \f$\ell_1,\ell_2,m\f$. This prefactor is defined by (cf Eq. (5.19))
 * \f[
 *      \Lambda_{\ell_1,\ell_2}^{(m)} = -\frac{2 N_{\ell_1,m} N_{\ell_2,m}}{\sqrt{\ell_1 (\ell_1+1) \ell_2 (\ell_2+1)}}
 * \f]
 *
 * If sign is not NULL, -1 will be stored in sign.
 *
 * The values are computed using the lgamma function to avoid overflows.
 *
 * Restrictions: \f$\ell_1,\ell_2 \ge 1\f$, \f$\ell_1,\ell_2 \ge m > 0\f$
 *
 * Symmetries: \f$\Lambda_{\ell_1,\ell_2}^{(m)} = \Lambda_{\ell_2,\ell_1}^{(m)}\f$
 *
 * This function is thread-safe.
 *
 * @param [in]  l1 \f$\ell_1\f$
 * @param [in]  l2 \f$\ell_2\f$
 * @param [in]  m  \f$m\f$
 * @param [out] sign set to -1 if sign != NULL
 * @retval lnLambda \f$\log{\Lambda_{\ell_1,\ell_2}^{(m)}}\f$
 */
edouble inline casimir_lnLambda(int l1, int l2, int m, sign_t *sign)
{
    if(sign != NULL)
        *sign = -1;
    return LOG2 + (loge(2.*l1+1)+loge(2*l2+1)-LOG4-loge(l1)-loge(l1+1)-loge(l2)-loge(l2+1)+lnfac(l1-m)+lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m))/2.0L;
}


/**
 * @brief Calculate \f$\epsilon(i\xi)\f$ for Drude model
 *
 * This function returns the dielectric function
 * \f[
 *      \epsilon(i\xi) = 1 + \frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in]  xi     \f$\xi\f$ imaginary frequency (in scaled units: \f$\xi=nT\f$)
 * @param [in]  omegap \f$\omega_\mathrm{P}\f$ Plasma frequency
 * @param [in]  gamma_ \f$\gamma\f$ relaxation frequency
 * @retval epsilon \f$\epsilon(\xi, \omega_\mathrm{P}, \gamma)\f$
 */
double casimir_epsilon(double xi, double omegap, double gamma_)
{
    return 1+omegap*omegap/(xi*(xi+gamma_));
}


/**
 * @brief Calculate \f$\log \epsilon(i\xi)\f$ for Drude model
 *
 * This function returns the logarithm of the dielectric function
 * \f[
 *      \epsilon(i\xi) = 1 + \frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in]  xi     \f$\xi\f$ imaginary frequency (in scaled units: \f$\xi=nT\f$)
 * @param [in]  omegap \f$\omega_\mathrm{P}\f$ Plasma frequency
 * @param [in]  gamma_ \f$\gamma\f$ relaxation frequency
 * @retval lnepsilon \f$\log{\epsilon(\xi, \omega_\mathrm{P}, \gamma)}\f$
 */
double casimir_lnepsilon(double xi, double omegap, double gamma_)
{
    return log1p(omegap*omegap/(xi*(xi+gamma_)));
}


/**
 * @brief Calculate Fresnel coefficients \f$r_{TE}\f$ and \f$r_{TM}\f$ for Drude model
 *
 * This function calculates the Fresnel coefficients for TE and TM mode
 *
 * This function is thread-safe.
 *
 * @param [in]      self    Casimir object
 * @param [in]      nT      \f$\xi=nT\f$ imaginary frequency
 * @param [in]      k       xy projection of wavevector
 * @param [in,out]  r_TE    Fresnel coefficient for TE mode
 * @param [in,out]  r_TM    Fresnel coefficient for TM mode
 */
void casimir_rp(casimir_t *self, edouble nT, edouble k, edouble *r_TE, edouble *r_TM)
{
    edouble epsilon = casimir_epsilon(nT, self->omegap_plane, self->gamma_plane);
    edouble beta = sqrte(1 + (epsilon-1)/(1 + pow_2(k/nT)));

    *r_TE = (1-beta)/(1+beta);
    *r_TM = (epsilon-beta)/(epsilon+beta);
}


/**
* @name converting
*/
/*@{*/

/**
 * @brief Convert free energy in SI units to free energy in units of \f$\mathcal{L}/(\hbar c)\f$
 *
 * This function returns 
 * \f[
 *      \mathcal{F}_\mathrm{scaled} = \mathcal{F}_\mathrm{SI} \frac{\mathcal{L}}{\hbar c}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in] F_SI free energy in SI units
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval free energy in scaled units
 */
double casimir_F_SI_to_scaled(double F_SI, double ScriptL)
{
    return ScriptL/(HBARC)*F_SI;
}


/**
 * @brief Convert free energy in units of \f$\hbar c / \mathcal{L}\f$ to SI units
 *
 * This function returns 
 * \f[
 *      \mathcal{F}_\mathrm{SI} = \mathcal{F}_\mathrm{scaled} \frac{\hbar c}{\mathcal{L}}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in] F_scaled \f$\mathcal{F}_\mathrm{scaled}\f$, free energy in units of \f$\hbar c / \mathcal{L}\f$
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval F_SI \f$\mathcal{F}_\mathrm{SI}\f$ free energy in units of Joule
 */
double casimir_F_scaled_to_SI(double F_scaled, double ScriptL)
{
    return HBARC/ScriptL*F_scaled;
}


/**
 * @brief Convert temperature in units of Kelvin to temperature in units of \f$\hbar c /(2\pi k_B \mathcal{L})\f$
 *
 * This function returns 
 * \f[
 *      T_\mathrm{scaled} = \frac{2\pi k_b \mathcal{L}}{\hbar c} T_\mathrm{SI}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in] T_SI temperature in units of Kelvin
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ (in units of meters)
 * @retval T_scaled temperature in units of \f$\hbar c /(2\pi k_B \mathcal{L})\f$
 */
double casimir_T_SI_to_scaled(double T_SI, double ScriptL)
{
    return 2*PI*KB*ScriptL/HBARC*T_SI;
}


/**
 * @brief Convert temperature in units of \f$\hbar c /(2\pi k_B \mathcal{L})\f$ to temperature in units of Kelvin
 *
 * This function returns 
 * \f[
 *      T_\mathrm{scaled} = \frac{2\pi k_b \mathcal{L}}{\hbar c} T_\mathrm{SI}
 * \f]
 *
 * This function is thread-safe.
 *
 * @param [in] T_scaled temperature in units of \f$\hbar c /(2\pi k_B \mathcal{L})\f$
 * @param [in] ScriptL \f$\mathcal{L} = R+L\f$ in units of meters
 * @retval T_SI temperature in units of Kelvin
 */
double casimir_T_scaled_to_SI(double T_scaled, double ScriptL)
{
    return HBARC/(2*PI*KB*ScriptL)*T_scaled;
}

/*@}*/


/**
 * @brief Calculate logarithm and sign of prefactor \f$\Xi_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of the prefactor for given
 * \f$\ell_1,\ell_2,m\f$. The prefactor is defined by Eq. (5.54).
 *
 * If sign is not NULL, the sign of \f$\Xi_{\ell_1 \ell_2}^{(m)}\f$ is stored in
 * sign.
 *
 * The values are computed using the lgamma function to avoid overflows.
 *
 * Restrictions: \f$\ell_1,\ell_2 \ge 1\f$, \f$\ell_1,\ell_2 \ge m > 0\f$
 *
 * This function is thread-safe.
 *
 * @param [in]  l1 \f$\ell_1\f$
 * @param [in]  l2 \f$\ell_2\f$
 * @param [in]  m  \f$m\f$
 * @param [out] sign
 * @retval logXi \f$\log{\Xi(\ell_1,\ell_2,m)}\f$
 */
edouble casimir_lnXi(int l1, int l2, int m, sign_t *sign)
{
    if(sign != NULL)
        *sign = MPOW(l2);
    return (loge(2*l1+1)+loge(2*l2+1)-lnfac(l1-m)-lnfac(l2-m)-lnfac(l1+m)-lnfac(l2+m)-loge(l1)-loge(l1+1)-loge(l2)-loge(l2+1))/2.0L \
           +lnfac(2*l1)+lnfac(2*l2)+lnfac(l1+l2)-LOG4*(2*l1+l2+1)-lnfac(l1-1)-lnfac(l2-1);
}


/**
* @name initialization and setting parameters
*/
/*@{*/

/**
 * @brief Create a new Casimir object for perfect reflectors
 *
 * This function will initialize a Casimir object with sphere and plane perfect
 * reflectors.
 *
 * Restrictions: \f$T > 0\f$, \f$0 < R/\mathcal{L} < 1\f$
 *
 * This function is not thread-safe.
 *
 * @param [out] self Casimir object
 * @param [in]  RbyScriptL \f$\frac{R}{\mathcal{L}} = \frac{R}{R+L}\f$
 * @param [in]  T temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 * @retval 0 if successful
 * @retval -1 if wrong value for RbyScriptL
 * @retval -2 if wrong value for T
 */
int casimir_init(casimir_t *self, double RbyScriptL, double T)
{
    double LbyR = 1./RbyScriptL - 1;
    if(RbyScriptL < 0 || RbyScriptL >= 1)
        return -1;
    if(T < 0)
        return -2;
    
    self->lmax = (int)ceil(CASIMIR_FACTOR_LMAX/LbyR);

    self->T           = T;
    self->RbyScriptL  = RbyScriptL;
    self->LbyR        = LbyR;
    self->precision   = CASIMIR_DEFAULT_PRECISION;
    self->verbose     = 0;
    self->cores       = 1;
    self->threads     = NULL;

    /* initialize mie cache */
    casimir_mie_cache_init(self);
    
    /* perfect reflectors */
    self->integration = -1; /* perfect reflectors */
    self->omegap_sphere = INFINITY;
    self->gamma_sphere  = 0;
    self->omegap_plane  = INFINITY;
    self->gamma_plane   = 0;

    return 0;
}


/**
 * @brief Set order of integration
 *
 * Set order/type of integration.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] integration: 0 perfect reflectors, >0: order of Gauss-Laguerre integration
 */
void casimir_set_integration(casimir_t *self, int integration)
{
    if(integration <= 0)
        self->integration = 0;
    else
        self->integration = integration;
}

/**
 * @brief Get order of integration
 *
 * Get order/type of integration.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval 0 if analytic integration for perfect reflectors
 * @retval order of integration for Gauss-Laguerre
 */
int casimir_get_integration(casimir_t *self)
{
    return self->integration;
}


/**
 * @brief Set \f$\omega_\mathrm{P}\f$ for the sphere
 *
 * Set the plasma frequency for the sphere.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] omegap \f$\omega_\mathrm{P}\f$ plasma frequency
 * @retval 1 if successful
 * @retval 0 if \f$\omega_\mathrm{P} < 0\f$
 */
int casimir_set_omegap_sphere(casimir_t *self, double omegap)
{
    if(omegap > 0)
    {
        self->omegap_sphere = omegap;
        self->integration   = 50;
        return 1;
    }
    return 0;
}

/**
 * @brief Set \f$\omega_\mathrm{P}\f$ for the plane
 *
 * Set the plasma frequency for the plane.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] omegap \f$\omega_\mathrm{P}\f$ plasma frequency
 * @retval 1 if successful
 * @retval 0 if \f$omega_\mathrm{P} < 0\f$
 */
int casimir_set_omegap_plane(casimir_t *self, double omegap)
{
    if(omegap > 0)
    {
        self->omegap_plane = omegap;
        self->integration  = 50;
        return 1;
    }
    return 0;
}


/**
 * @brief Get \f$\omega_\mathrm{P}\f$ for the sphere
 *
 * Get the plasma frequency for the sphere.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval plasma frequency
 */
double casimir_get_omegap_sphere(casimir_t *self)
{
    return self->omegap_sphere;
}

/**
 * @brief Get \f$\omega_\mathrm{P}\f$ for the plane
 *
 * Get the plasma frequency \f$\omega_\mathrm{P}\f$ for the plane.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval omegap \f$\omega_\mathrm{P}\f$plasma frequency
 */
double casimir_get_omegap_plane(casimir_t *self)
{
    return self->omegap_plane;
}


/**
 * @brief Set \f$\gamma\f$ for the sphere
 *
 * Set the relaxation frequency for the sphere.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] gamma_ \f$\gamma\f$ relaxation frequency
 * @retval 1 if successful
 * @retval 0 if \f$\gamma < 0\f$
 */
int casimir_set_gamma_sphere(casimir_t *self, double gamma_)
{
    if(gamma_ > 0)
    {
        self->gamma_sphere = gamma_;
        self->integration = 50;
        return 1;
    }
    return 0;
}

/**
 * @brief Set \f$\gamma\f$ for the plane
 *
 * Set the relaxation frequency \f$\gamma\f$ for the plane.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] gamma_ relaxation frequency
 * @retval 1 if successful
 * @retval 0 if \f$\gamma < 0\f$
 */
int casimir_set_gamma_plane(casimir_t *self, double gamma_)
{
    if(gamma_ > 0)
    {
        self->gamma_plane = gamma_;
        self->integration = 50;
        return 1;
    }
    return 0;
}


/**
 * @brief Get \f$\gamma\f$ for the sphere
 *
 * Get the relaxation frequency \f$\gamma\f$ for the sphere.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval gamma \f$\gamma\f$ relaxation frequency
 */
double casimir_get_gamma_sphere(casimir_t *self)
{
    return self->gamma_sphere;
}

/**
 * @brief Get \f$\gamma\f$ for the plane
 *
 * Get the relaxation frequency for the plane.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval gamma \f$\gamma\f$ relaxation frequency
 */
double casimir_get_gamma_plane(casimir_t *self)
{
    return self->gamma_plane;
}


/**
 * @brief Return numbers of used cores
 *
 * See casimir_set_cores.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval number of used cores (>=0)
 */
int casimir_get_cores(casimir_t *self)
{
    return self->cores;
}


/**
 * @brief Set the number of used cores
 *
 * This library supports multiple processor cores. However, you must specify
 * how many cores the library should use. By default, only one core will be
 * used.
 *
 * The libraray uses POSIX threads for parallelization. Threads share memory and
 * for this reason all cores must be on the same computer.
 *
 * Restrictions: cores > 0
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] cores number of cores that should be used
 * @retval 1 if successful
 * @retval 0 if cores < 1
 */
int casimir_set_cores(casimir_t *self, int cores)
{
    if(cores < 1)
        return 0;

    self->cores = cores;
    self->threads = xrealloc(self->threads, cores*sizeof(pthread_t));

    return 1;
}


/**
 * @brief Set maximum value of l
 *
 * In general the round trip matrices are infinite. For a numerical evaluation
 * the dimension has to be limited to a finite value. The accuracy of the
 * result depends on the truncation of the vector space. For more information,
 * cf. chapter 6.1.
 *
 * Please note that the cache of the Mie coefficients has to be cleaned. This
 * means that all Mie coefficients have to be calculated again.
 *
 * In order to get meaningful results and to prevent recalculating Mie
 * coefficients, set the lmax at the beginning before doing expensive
 * computations.
 *
 * This function is not thread-safe.
 * 
 * @param [in,out] self Casimir object
 * @param [in] lmax maximum number of \f$\ell\f$
 * @retval 1 if successful
 * @retval 0 if lmax < 1
 */
int casimir_set_lmax(casimir_t *self, int lmax)
{
    if(lmax <= 0)
        return 0;

    self->lmax = lmax;

    /* reinit mie cache */
    casimir_mie_cache_clean(self);

    return 1;
}


/**
 * @brief Get maximum value of l
 *
 * See casimir_set_lmax.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval lmax maximum value of \f$\ell\f$
 */
int casimir_get_lmax(casimir_t *self)
{
    return self->lmax;
}


/**
 * @brief Get verbose flag
 *
 * Return if the verbose flag is set.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval 0 if verbose flag is not set
 * @retval 1 if verbose flag is set
 */
int casimir_get_verbose(casimir_t *self)
{
    return self->verbose;
}


/**
 * @brief Set verbose flag
 *
 * Use this function to set the verbose flag. If set to 1, this will cause the
 * library to print information to stderr. 
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] verbose 1 if verbose, 0 if not verbose
 * @retval 1
 */
int casimir_set_verbose(casimir_t *self, int verbose)
{
    self->verbose = verbose ? 1 : 0;
    return 1;
}


/**
 * @brief Get precision
 *
 * See casimir_set_precision
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @retval precision \f$\epsilon_p\f$
 */
double casimir_get_precision(casimir_t *self)
{
    return self->precision;
}


/**
 * @brief Set precision
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] precision \f$\epsilon_p\f$
 * @retval 1 if successful
 * @retval 0 if precision <= 0
 */
int casimir_set_precision(casimir_t *self, double precision)
{
    if(precision <= 0)
        return 0;

    self->precision = precision;
    return 1;
}


/**
 * @brief Free memory for Casimir object
 *
 * This function will free allocated memory for the Casimir object. If you have
 * allocated memory for the object yourself, you have, however, to free this
 * yourself.
 *
 * This function is not thread-safe.
 * 
 * @param [in,out] self Casimir object
 */
void casimir_free(casimir_t *self)
{
    casimir_mie_cache_free(self);

    if(self->threads != NULL)
    {
        xfree(self->threads);
        self->threads = NULL;
    }
}

/*@}*/


/**
* @name Mie coefficients
*/
/*@{*/

/** Return the logarithm of the prefactors \f$a_{\ell,0}^\mathrm{perf}\f$, \f$b_{\ell,0}^\mathrm{perf}\f$ and their signs
 *
 * For small frequencies \f$\chi = \frac{\xi R}{c} \ll 1\f$ the Mie
 * coeffiecients scale like
 * \f[
 * a_{\ell}^\mathrm{perf} = a_{\ell,0}^\mathrm{perf} \left(\frac{\chi}{2}\right)^{2\ell+1} \\
 * \f]
 * \f[
 * b_{\ell}^\mathrm{perf} = b_{\ell,0}^\mathrm{perf} \left(\frac{\chi}{2}\right)^{2\ell+1}
 * \f]
 * This function returns the logarithm of the prefactors
 * \f$a_{\ell,0}^\mathrm{perf}\f$, \f$b_{\ell,0}^\mathrm{perf}\f$ and their
 * signs.
 *
 * In scaled units: \f$\chi = nT \frac{R}{\mathcal{L}}\f$
 *
 * This function is thread-safe.
 *
 * @param [in] order \$\ell\f$
 * @param [out] a0 coefficient \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_a0 sign of \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] b0 coefficient \f$b_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_b0 sign of \f$b_{\ell,0}^\mathrm{perf}\f$
 */
void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0)
{
    *sign_a0 = MPOW(l);
    *sign_b0 = MPOW(l+1);
    *b0 = LOGPI-lngamma(l+0.5)-lngamma(l+1.5);
    *a0 = *b0+log1p(1.0/l);
}


/**
 * @brief Return logarithm of Mie coefficient \f$a_\ell\f$ for perfect reflectors and its sign
 *
 * The frequency will be determined by n: \f$\xi = nT\f$
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * This function is thread-safe - as long you don't change temperature and
 * aspect ratio.
 *
 * @param [in,out] self Casimir object
 * @param [in] order \f$\ell\f$
 * @param [in] n Matsubara term, \f$xi = nT\f$
 * @param [out] sign sign of \f$a_\ell\f$
 * @retval logarithm of Mie coefficient \f$a_\ell\f$
 */
double casimir_lna_perf(casimir_t *self, const int l, const int n, sign_t *sign)
{
    edouble nominator, denominator, frac, ret;
    edouble lnKlp,lnKlm,lnIlm,lnIlp;
    edouble prefactor;
    edouble chi = n*self->T*self->RbyScriptL;
    edouble lnfrac = log(chi)-log(l);

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother? 
     */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp);

    prefactor = LOGPI-LOG2+lnIlp-lnKlp;
    *sign = MPOW(l+1);

    /* numinator */
    {
        frac = expe(lnfrac+lnIlm-lnIlp);
        if(frac < 1)
            nominator = log1pe(fabse(-frac));
        else
        {
            if(frac > 1)
                *sign *= -1;

            nominator = loge(fabse(1-frac));
        }
    }
    /* denominator */
    {
        frac = expe(lnfrac+lnKlm-lnKlp);
        if(frac < 1)
            denominator = log1pe(frac);
        else
            denominator = log1pe(frac);
    }

    ret = prefactor+nominator-denominator;

    return ret;
}


/**
 * @brief Return logarithm of Mie coefficient \f$b_\ell\f$ for perfect reflectors and its sign
 *
 * The frequency will be determined by n: \f$\xi = nT\f$
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * This function is thread safe - as long you don't change temperature and
 * aspect ratio.
 *
 * @param [in,out] self Casimir object
 * @param [in] order \f$\ell\f$
 * @param [in] n Matsubara term, \f$xi = nT\f$
 * @param [out] sign sign of \f$b_\ell\f$
 * @retval logarithm of Mie coefficient \f$b_\ell\f$
 */
double casimir_lnb_perf(casimir_t *self, const int l, const int n, sign_t *sign)
{
    edouble chi = n*self->T*self->RbyScriptL;
    edouble lnInu, lnKnu;

    bessel_lnInuKnu(l, chi, &lnInu, &lnKnu);
    *sign = MPOW(l+1);

    return LOGPI-LOG2+lnInu-lnKnu;
}


/**
 * @brief Return logarithm of Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for Drude model
 *
 * For \f$\omega_\mathrm{P} = \infty\f$ the Mie coefficient for perfect
 * reflectors are returned (see casimir_lna_perf and casimir_lnb_perf).
 *
 * Cf. Eqs. (3.30) and (3.31).
 *
 * sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * This function is thread safe - as long you don't change temperature, aspect
 * ratio and dielectric properties of sphere.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term, \f$\xi = nT\f$
 * @param [in] order \f$\ell\f$
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 * @param [out] sign_a sign of Mie coefficient \f$a_\ell\f$
 * @param [out] sign_b sign of Mie coefficient \f$b_\ell\f$
 */
void casimir_lnab(casimir_t *self, const int n_mat, const int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
{ 
    sign_t sign_sla, sign_slb, sign_slc, sign_sld;
    edouble ln_n, ln_sla, ln_slb, ln_slc, ln_sld;
    edouble lnIl, lnKl, lnIlm, lnKlm, lnIl_nchi, lnKl_nchi, lnIlm_nchi, lnKlm_nchi;
    edouble xi = n_mat*self->T;
    edouble chi = xi*self->RbyScriptL;
    edouble ln_chi = loge(xi)+loge(self->RbyScriptL);
    edouble omegap = self->omegap_sphere;
    edouble gamma_ = self->gamma_sphere;
    sign_t sign_a_num, sign_a_denom, sign_b_num, sign_b_denom;

    if(isinf(omegap))
    {
        *lna = casimir_lna_perf(self, l, n_mat, sign_a);
        *lnb = casimir_lnb_perf(self, l, n_mat, sign_b);
        return;
    }

    ln_n = casimir_lnepsilon(xi, omegap, gamma_)/2;

    bessel_lnInuKnu(l,   chi, &lnIl,  &lnKl);
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);

    bessel_lnInuKnu(l,   expe(ln_n)*chi, &lnIl_nchi,  &lnKl_nchi);
    bessel_lnInuKnu(l-1, expe(ln_n)*chi, &lnIlm_nchi, &lnKlm_nchi);

    ln_sla = lnIl_nchi + logadd_s(lnIl,      +1, ln_chi+lnIlm,           -1, &sign_sla);
    ln_slb = lnIl      + logadd_s(lnIl_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_slb);
    ln_slc = lnIl_nchi + logadd_s(lnKl,      +1, ln_chi+lnKlm,           +1, &sign_slc);
    ln_sld = lnKl      + logadd_s(lnIl_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_sld);

    /* XXX FIXME XXX */
    /*
    printf("n =%.15g\n", (double)expe(ln_n));
    printf("n2=%.15g\n", (double)expe(2*ln_n));
    printf("lnIl = %.15g\n", (double)expe(lnIl));
    printf("chi=%.15g\n", (double)chi);
    */

    /*
    printf("sla=%.15g\n", (double)(sign_sla*expe(ln_sla)));
    printf("slb=%.15g\n", (double)(sign_slb*expe(ln_slb)));
    printf("slc=%.15g\n", (double)(sign_slc*expe(ln_slc)));
    printf("sld=%.15g\n", (double)(sign_sld*expe(ln_sld)));
    */

    *lna = LOGPI - LOG2 + logadd_s(2*ln_n+ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_a_num) - logadd_s(2*ln_n+ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_a_denom);
    *lnb = LOGPI - LOG2 + logadd_s(       ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_b_num) - logadd_s(       ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_b_denom);

    *sign_a = sign_a_num*sign_a_denom;
    *sign_b = sign_b_num*sign_b_denom;
}


/**
 * @brief Initialize Mie cache
 *
 * This will initialize a cache for the Mie coefficients. If a Mie coefficient
 * \f$a_\ell(\xi=nT)\f$ or \f$b_\ell(xi=nT)\f$ is needed, the coefficients only
 * needs to be calculated the first time. The result will be saved.
 *
 * The cache will grow on demand. But it will never shrink - unless you call
 * casimir_free or casimir_mie_cache_clean.
 *
 * The casimir_init function will call this function and the cache will be used
 * throughout computations. So you usually don't want to use this function
 * yourself.
 *
 * This function is not thread safe.
 *
 * @param [in,out] casimir object
 */
void casimir_mie_cache_init(casimir_t *self)
{
    casimir_mie_cache_t *cache = self->mie_cache = (casimir_mie_cache_t *)xmalloc(sizeof(casimir_mie_cache_t));

    cache->lmax = self->lmax;
    cache->nmax = 0;

    pthread_mutex_init(&cache->mutex, NULL);

    cache->entries = xmalloc(sizeof(casimir_mie_cache_entry_t *));

    cache->entries[0] = xmalloc(sizeof(casimir_mie_cache_entry_t));
    cache->entries[0]->ln_al   = NULL;
    cache->entries[0]->sign_al = NULL;
    cache->entries[0]->ln_bl   = NULL;
    cache->entries[0]->sign_bl = NULL;
}


/**
 * @brief Allocate memory for the Mie coefficients \f$a_\ell\f$ and \f$b_\ell\f$
 *
 * This function computes all the Mie coefficients for \f$\xi=nT\f$ and stores
 * it in cache. Make sure you have already initialized the cache (cf.
 * casimir_mie_cache_init). Initialization is done by casimir_init.
 *
 * This function will be called by casimir_mie_cache_get if a Mie coefficient
 * is not precomputed yet.
 *
 * You usually don't want to use this function yourself.
 *
 * This function is not thread safe.
 *
 * @param [in,out] self Casimir object
 * @param [in,out] cache Mie cache
 */
void casimir_mie_cache_alloc(casimir_t *self, int n)
{
    int l;
    const int nmax = self->mie_cache->nmax;
    const int lmax = self->mie_cache->lmax;
    casimir_mie_cache_t *cache = self->mie_cache;

    if(n > nmax)
    {
        cache->entries = xrealloc(cache->entries, (n+1)*sizeof(casimir_mie_cache_entry_t *));

        for(l = nmax+1; l <= n; l++)
            cache->entries[l] = NULL;
    }

    if(cache->entries[n] == NULL)
    {
        cache->entries[n] = xmalloc(sizeof(casimir_mie_cache_entry_t));
        casimir_mie_cache_entry_t *entry = cache->entries[n];

        entry->ln_al   = xmalloc((lmax+1)*sizeof(double));
        entry->ln_bl   = xmalloc((lmax+1)*sizeof(double));
        entry->sign_al = xmalloc((lmax+1)*sizeof(sign_t));
        entry->sign_bl = xmalloc((lmax+1)*sizeof(sign_t));

        entry->ln_al[0] = entry->ln_bl[0] = NAN; /* should never be read */
        for(l = 1; l <= lmax; l++)
            casimir_lnab(self, n, l, &entry->ln_al[l], &entry->ln_bl[l], &entry->sign_al[l], &entry->sign_bl[l]);
    }

    self->mie_cache->nmax = n;
}


/**
 * @brief Clean memory of cache
 *
 * This function will free allocated memory for the cache. The cache will still
 * work, but the precomputed values for the Mie coefficients will be lost.
 *
 * You usually don't want to use this function yourself.
 *
 * This function is not thread-safe.
 *
 * @param [in, out] Casimir object
 */
void casimir_mie_cache_clean(casimir_t *self)
{
    casimir_mie_cache_free(self);
    casimir_mie_cache_init(self);
}

/**
 * @brief Clean memory of cache
 *
 * Get Mie coefficients for \f$\ell\f$ and Matsubara frequency \f$\xi=nT\f$. If
 * the Mie coefficients have not been calculated yet, they will be computed and
 * stored in the cache.
 *
 * This function is thread-safe.
 *
 * @param [in, out] Casimir object
 * @param [in] order \f$\ell\f$
 * @param [in] n Matsubara term
 * @param [out] ln_a logarithm of \f$a_\ell\f$
 * @param [out] sign_a sign of \f$a_\ell\f$
 * @param [out] ln_b logarithm of \f$b_\ell\f$
 * @param [out] sign_b sign of \f$b_\ell\f$
 */
void casimir_mie_cache_get(casimir_t *self, int l, int n, double *ln_a, sign_t *sign_a, double *ln_b, sign_t *sign_b)
{
    casimir_mie_cache_entry_t *entry;
    int nmax;

    /* this mutex is important to prevent memory corruption */
    pthread_mutex_lock(&self->mie_cache->mutex);

    nmax = self->mie_cache->nmax;
    if(n > nmax || self->mie_cache->entries[n] == NULL)
        casimir_mie_cache_alloc(self, n);

    /* This is a first class example of concurrent accesses on memory and
     * locking: You might think, we don't need a lock anymore. All the data has
     * been written and it is safe to release the mutex. Unfortunately, this is
     * not true.
     *
     * Although all data has been written, another thread might add more Mie
     * coefficients to the cache. In this case, the array
     * self->mie_cache->entries needs to hold more values. For this reason,
     * realloc is called. However, realloc may change the position of the array
     * in memory and the pointer self->mie_cache->entries becomes invalid. This
     * will usually cause a segmentation fault.
     */
    entry = self->mie_cache->entries[n];

    /* at this point it is finally is safe to release the mutex without lock */
    pthread_mutex_unlock(&self->mie_cache->mutex);

    *ln_a   = entry->ln_al[l];
    *sign_a = entry->sign_al[l];
    *ln_b   = entry->ln_bl[l];
    *sign_b = entry->sign_bl[l];
}

/**
 * @brief Free memory of cache.
 *
 * This function will free allocated memory for the cache. Don't use the cache
 * afterwards!
 *
 * You usually don't want to use this function yourself.
 *
 * This function is not thread-safe.
 *
 * @param [in, out] Casimir object
 */
void casimir_mie_cache_free(casimir_t *self)
{
    int n;
    casimir_mie_cache_t *cache = self->mie_cache;
    casimir_mie_cache_entry_t **entries = cache->entries;

    pthread_mutex_destroy(&cache->mutex);

    /* free
     * 1) the lists of al, bl, sign_al, sign_bl for every entry
     * 2) every entry (casimir_mie_cache_entry_t)
     * 3) the list of entries
     * 4) the mie cache object (casimir_mie_cache_t)
     */
    for(n = 0; n <= cache->nmax; n++)
    {
        if(entries[n] != NULL)
        {
            if(entries[n]->ln_al != NULL)
                xfree(entries[n]->ln_al);
            if(entries[n]->sign_al != NULL)
                xfree(entries[n]->sign_al);
            if(entries[n]->ln_bl != NULL)
                xfree(entries[n]->ln_bl);
            if(entries[n]->sign_bl != NULL)
                xfree(entries[n]->sign_bl);

            xfree(entries[n]);
        }
    }

    xfree(cache->entries);
    cache->entries = NULL;

    xfree(self->mie_cache);
    self->mie_cache = NULL;
}

/*@}*/

/* Sum len numbers in value.
   The idea is: To avoid a loss of significance, we sum beginning with smallest
   number and add up in increasing order
*/
static double _sum(double values[], size_t len)
{
    int i;
    double sum = 0;

    for(i = len-1; i > 0; i--)
        sum += values[i];

    sum += values[0]/2;

    return sum;
}

/* This is the function the thread executes */
static void *_thread_f(void *p)
{
    casimir_thread_t *r = (casimir_thread_t *)p;
    r->value = casimir_F_n(r->self, r->n, &r->nmax);
    return r;
}

/* This function starts a thread to compute the free energy corresponding to
 * the Matsubara term n
 *
 * The memory for the thread object and the variable r will be allocated in this function
 * and freed in _join_threads.
 */
static pthread_t *_start_thread(casimir_t *self, int n)
{
    pthread_t *t = xmalloc(sizeof(pthread_t));
    casimir_thread_t *r = xmalloc(sizeof(casimir_thread_t));

    r->self  = self;
    r->n     = n;
    r->value = 0;
    r->nmax  = 0;

    pthread_create(t, NULL, _thread_f, (void *)r);

    return t;
}

/* This function tries to join threads and writed the result to values */
static int _join_threads(casimir_t *self, double values[], int *ncalc)
{
    int i, joined = 0, running = 0;
    casimir_thread_t *r;
    pthread_t **threads = self->threads;

    for(i = 0; i < self->cores; i++)
    {
        if(threads[i] != NULL)
        {
            running++;

            if(pthread_tryjoin_np(*threads[i], (void *)&r) == 0)
            {
                joined++;

                if(r->n > *ncalc)
                    *ncalc = r->n;

                values[r->n] = r->value;
                xfree(r);
                xfree(threads[i]);
                threads[i] = NULL;
            }
        }
    }

    if(running == 0)
        return -1;

    return joined;
}


/**
* @name Calculate free energy
*/
/*@{*/

/**
 * @brief Calculate free energy for Matsubara term n
 *
 * This function calculates the free energy for the Matsubara term n. If mmax
 * is not NULL, the largest value of m calculated will be stored in mmax.
 *
 * This function is thread-safe - as long you don't change temperature, aspect
 * ratio, dielectric properties of sphere and plane, lmax, integration and
 * verbose.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term, \f$\xi=nT\f$
 * @param [out] mmax maximum number of \f$m\f$
 * @retval F Casimir free energy for given \f$n\f$
 */
double casimir_F_n(casimir_t *self, const int n, int *mmax)
{
    double precision = self->precision;
    double sum_n = 0;
    int m;
    const int lmax = self->lmax;
    double values[lmax+1];
    integration_perf_t int_perf;

    /* perfect reflectors */

    for(m = 0; m <= lmax; m++)
        values[m] = 0;

    if(self->integration <= 0) 
        casimir_integrate_perf_init(&int_perf, n*self->T, self->lmax);

    for(m = 0; m <= self->lmax; m++)
    {
        values[m] = casimir_logdetD(self,n,m,&int_perf);

        if(self->verbose)
            fprintf(stderr, "# n=%d, m=%d, value=%.15g\n", n, m, values[m]);

        /* If F is !=0 and value/F < 1e-16, then F+value = F. The addition
         * has no effect.
         * As for larger m value will be even smaller, we can skip the
         * summation here. 
         */
        sum_n = _sum(values, lmax+1);
        if(values[0] != 0 && fabs(values[m]/sum_n) < precision)
            break;
    }

    if(self->integration <= 0) 
        casimir_integrate_perf_free(&int_perf);

    if(self->verbose)
        fprintf(stderr, "# n=%d, value=%.15g\n", n, sum_n);

    if(mmax != NULL)
        *mmax = m;

    return sum_n;
}


/**
 * @brief Calculate free energy
 *
 * This function calculates the free energy. If nmax is not NULL, the highest
 * Matsubara term calculated will be stored in nnmax.
 *
 * This function will use as many cores as specified by casimir_set_cores.
 * 
 * @param [in,out] self Casimir object
 * @param [out] nmax maximum number of n
 * @retval F Casimir free energy
 */
double casimir_F(casimir_t *self, int *nmax)
{
    int i, n = 0;
    double sum_n = 0;
    const double precision = self->precision;
    double *values = NULL;
    size_t len = 0;
    int ncalc = 0;
    const int cores = self->cores;
    const int delta = MAX(1024, cores);
    pthread_t **threads = self->threads;

    for(i = 0; i < cores; i++)
        threads[i] = NULL;

    /* So, here we sum up all m and n that contribute to F.
     * So, what do we do here?
     *
     * We want to evaluate
     *      \sum_{n=0}^\infty \sum{m=0}^{l_max} log(det(1-M)),
     * where the terms for n=0 and m=0 are weighted with a factor 1/2.
     */
    while(1)
    {
        if(n >= len)
        {
            values = (double *)xrealloc(values, (len+delta)*sizeof(double));

            for(i = len; i < len+delta; i++)
                values[i] = 0;

            len += delta;
        }

        if(cores > 1)
        {
            for(i = 0; i < cores; i++)
                if(threads[i] == NULL)
                    threads[i] = _start_thread(self, n++);

            if(_join_threads(self, values, &ncalc) == 0)
                usleep(CASIMIR_IDLE);
        }
        else
        {
            values[n] = casimir_F_n(self, n, NULL);

            ncalc = n;
            n++;
        }

        if(values[0] != 0)
        {
            if(fabs(values[ncalc]/(2*values[0])) < precision)
            {
                if(cores > 1)
                    while(_join_threads(self, values, &ncalc) != -1)
                        usleep(CASIMIR_IDLE);

                sum_n = _sum(values, n);

                /* get out of here */
                if(nmax != NULL)
                    *nmax = n-1; // we calculated n term from n=0,...,nmax=n-1

                if(values != NULL)
                    xfree(values);

                return self->T/PI*sum_n;
            }
        }
    }
}


/**
 * @brief Calculate \f$\log\det \mathcal{D}^{(m)}(\xi=0)\f$ for EE and MM
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n=0\f$.
 *
 * This function is thread-safe as long you don't change lmax and aspect ratio.
 *
 * @param [in,out] self Casimir object
 * @param [in] m
 * @param [out] logdet_EE
 * @param [out] logdet_MM
 */
void casimir_logdetD0(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    int l1,l2,min,max,dim;
    const edouble lnRbyScriptL = loge(self->RbyScriptL);
    matrix_edouble_t *EE = NULL, *MM = NULL;
    matrix_sign_t *EE_sign = NULL, *MM_sign = NULL;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    if(logdet_EE != NULL)
    {
        EE = matrix_edouble_alloc(dim);
        EE_sign = matrix_sign_alloc(dim);
    }
    if(logdet_MM != NULL)
    {
        MM = matrix_edouble_alloc(dim);
        MM_sign = matrix_sign_alloc(dim);
    }

    /* calculate the logarithm of the matrix elements of D */
    for(l1 = min; l1 <= max; l1++)
        for(l2 = min; l2 <= max; l2++)
        {
            /* i: row of matrix, j: column of matrix */
            const int i = l1-min, j = l2-min;
            sign_t sign_a0, sign_b0, sign_xi;
            double lna0, lnb0;
            const edouble lnXiRL = casimir_lnXi(l1,l2,m,&sign_xi)+(2*l1+1)*lnRbyScriptL;
            casimir_lnab0(l1, &lna0, &sign_a0, &lnb0, &sign_b0);
            edouble v;
            sign_t sign;

            if(EE != NULL)
            {
                if(l1 != l2)
                {
                    v    = lna0+lnXiRL;
                    sign = -sign_xi*sign_a0;
                }
                else
                    v = logadd_s(0, +1, lna0+lnXiRL, -sign_xi*sign_a0, &sign);

                matrix_set(EE, i,j, v);
                matrix_set(EE_sign, i,j, sign);

                TERMINATE(isinf(v), "EE l1=%d,l2=%d: inf (lna0=%g, lnXiRL=%g)", l1, l2, lna0, (double)lnXiRL);
                TERMINATE(isnan(v), "EE l1=%d,l2=%d: nan (lna0=%g, lnXiRL=%g)", l1, l2, lna0, (double)lnXiRL);
            }
            if(MM != NULL)
            {
                if(l1 != l2)
                {
                    v    = lnb0+lnXiRL;
                    sign = sign_xi*sign_b0;
                }
                else
                    v = logadd_s(0, +1, lnb0+lnXiRL, sign_xi*sign_b0, &sign);

                matrix_set(MM, i,j, v);
                matrix_set(MM_sign, i,j, sign);

                TERMINATE(isinf(v), "MM l1=%d,l2=%d: inf (lnb0=%g, lnXiRL=%g)", l1, l2, lnb0, (double)lnXiRL);
                TERMINATE(isnan(v), "MM l1=%d,l2=%d: nan (lnb0=%g, lnXiRL=%g)", l1, l2, lnb0, (double)lnXiRL);
            }
        }

    /* balance matrices, calculate logdet and free space */
    if(EE != NULL)
    {
        matrix_edouble_log_balance(EE);

        #ifdef USE_LAPACK
            *logdet_EE = matrix_logdet_lapack(EE, EE_sign);
        #else
            matrix_edouble_exp(EE, EE_sign);
            *logdet_EE = matrix_edouble_logdet(EE);
        #endif

        matrix_sign_free(EE_sign);
        matrix_edouble_free(EE);
    }
    if(MM != NULL)
    {
        matrix_edouble_log_balance(MM);

        #ifdef USE_LAPACK
            *logdet_MM = matrix_logdet_lapack(MM, MM_sign);
        #else
            matrix_edouble_exp(MM, MM_sign);
            *logdet_MM = matrix_edouble_logdet(MM);
        #endif

        matrix_sign_free(MM_sign);
        matrix_edouble_free(MM);
    }
}


/**
 * @brief Calculate \f$\log\det \mathcal{D}^{(m)}(\xi=nT)\f$
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n\f$.
 *
 * This function is thread-safe - as long you don't change lmax, temperature,
 * aspect ratio, dielectric properties of sphere and plane, and integration.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term
 * @param [in] m
 * @param [in,out] m integration object
 * @retval logdetD \f$\log \det \mathcal{D}^{(m)}(\xi=nT)\f$
 */
double casimir_logdetD(casimir_t *self, int n, int m, void *integration_obj)
{
    int min,max,dim,l1,l2;
    double logdet = 0;
    double nTRbyScriptL = n*self->T*self->RbyScriptL;

    min = MAX(m,1);
    max = self->lmax;

    dim = (max-min+1);

    if(n == 0)
    {
        double logdet_EE = 0, logdet_MM = 0;

        if(isinf(self->omegap_plane))
            /* perfect reflector */
            casimir_logdetD0(self, m, &logdet_EE, &logdet_MM);
        else
            /* drude mirror */
            casimir_logdetD0(self, m, &logdet_EE, NULL);

        return logdet_EE + logdet_MM;
    }

    matrix_edouble_t *M = matrix_edouble_alloc(2*dim);
    matrix_sign_t *M_sign = matrix_sign_alloc(2*dim);

    /* M_EE, -M_EM
       M_ME,  M_MM */
    for(l1 = min; l1 <= max; l1++)
    {
        for(l2 = min; l2 <= l1; l2++)
        {
            const int Delta_ij = (l1 == l2 ? 0 : -INFINITY);
            const int i = l1-min, j = l2-min;
            casimir_integrals_t cint;
            double ln_al1, ln_bl1, ln_al2, ln_bl2;
            sign_t sign_al1, sign_bl1, sign_al2, sign_bl2;

            casimir_mie_cache_get(self, l1, n, &ln_al1, &sign_al1, &ln_bl1, &sign_bl1);
            casimir_mie_cache_get(self, l2, n, &ln_al2, &sign_al2, &ln_bl2, &sign_bl2);

            if(nTRbyScriptL < 1)
            {
                double lognTRbyScriptL = log(nTRbyScriptL);
                ln_al1 -= (l1-l2)*lognTRbyScriptL;
                ln_bl1 -= (l1-l2)*lognTRbyScriptL;

                ln_al2 -= (l2-l1)*lognTRbyScriptL;
                ln_bl2 -= (l2-l1)*lognTRbyScriptL;
            }

            if(self->integration > 0)
                casimir_integrate_drude(self, &cint, l1, l2, m, n, self->T);
            else
                casimir_integrate_perf(integration_obj, l1, l2, m, &cint);

            /* EE */
            {
                sign_t sign;
                edouble list_ij[] = { Delta_ij, ln_al1+cint.lnA_TE, ln_al1+cint.lnB_TM };
                edouble list_ji[] = { Delta_ij, ln_al2+cint.lnA_TE, ln_al2+cint.lnB_TM };

                sign_t signs_ij[] = { +1, -            sign_al1*cint.signA_TE, -            sign_al1*cint.signB_TM };
                sign_t signs_ji[] = { +1, -MPOW(l1+l2)*sign_al2*cint.signA_TE, -MPOW(l1+l2)*sign_al2*cint.signB_TM };

                matrix_set(M, i,j, logadd_ms(list_ij, signs_ij, 3, &sign));
                matrix_set(M_sign, i,j, sign);

                matrix_set(M, j,i, logadd_ms(list_ji, signs_ji, 3, &sign));
                matrix_set(M_sign, j,i, sign);

                TERMINATE(isinf(matrix_get(M,i,j)), "EE, i=%d,j=%d is inf", i,j);
                TERMINATE(isnan(matrix_get(M,i,j)), "EE, i=%d,j=%d is nan", i,j);

                TERMINATE(isinf(matrix_get(M,j,i)), "EE, i=%d,j=%d is inf", j,i);
                TERMINATE(isnan(matrix_get(M,j,i)), "EE, i=%d,j=%d is nan", j,i);
            }

            /* MM */
            {
                sign_t sign;
                edouble list_ij[] = { Delta_ij, ln_bl1+cint.lnA_TM, ln_bl1+cint.lnB_TE };
                edouble list_ji[] = { Delta_ij, ln_bl2+cint.lnA_TM, ln_bl2+cint.lnB_TE };

                sign_t signs_ij[] = { +1, -            sign_bl1*cint.signA_TM, -            sign_bl1*cint.signB_TE };
                sign_t signs_ji[] = { +1, -MPOW(l1+l2)*sign_bl2*cint.signA_TM, -MPOW(l1+l2)*sign_bl2*cint.signB_TE };

                matrix_set(M, i+dim,j+dim, logadd_ms(list_ij, signs_ij, 3, &sign));
                matrix_set(M_sign, i+dim,j+dim, sign);

                matrix_set(M, j+dim,i+dim, logadd_ms(list_ji, signs_ji, 3, &sign));
                matrix_set(M_sign, j+dim,i+dim, sign);

                TERMINATE(isinf(matrix_get(M,i+dim,j+dim)), "MM, i=%d,j=%d is inf", i+dim,j+dim);
                TERMINATE(isnan(matrix_get(M,i+dim,j+dim)), "MM, i=%d,j=%d is nan", i+dim,j+dim);

                TERMINATE(isinf(matrix_get(M,j+dim,i+dim)), "MM, i=%d,j=%d is inf", j+dim,i+dim);
                TERMINATE(isnan(matrix_get(M,j+dim,i+dim)), "MM, i=%d,j=%d is nan", j+dim,i+dim);
            }


            if(m != 0)
            {
                /* M_EM */
                {
                    sign_t sign;
                    edouble list_ij[] = { ln_al1+cint.lnC_TE, ln_al1+cint.lnD_TM };
                    sign_t signs_ij[] = { -sign_al1*cint.signC_TE,  -sign_al1*cint.signD_TM };

                    edouble list_ji[] = { ln_al2+cint.lnD_TE, ln_al2+cint.lnC_TM };
                    sign_t signs_ji[] = { -MPOW(l1+l2+1)*sign_al2*cint.signD_TE, -MPOW(l1+l2+1)*sign_al2*cint.signC_TM };

                    matrix_set(M, dim+i,j, logadd_ms(list_ij, signs_ij, 2, &sign));
                    matrix_set(M_sign, dim+i,j, sign);

                    matrix_set(M, dim+j,i, logadd_ms(list_ji, signs_ji, 2, &sign));
                    matrix_set(M_sign, dim+j,i, sign);

                    TERMINATE(isinf(matrix_get(M,i+dim,j)), "EM, i=%d,j=%d is inf", i+dim,j);
                    TERMINATE(isnan(matrix_get(M,i+dim,j)), "EM, i=%d,j=%d is nan", i+dim,j);

                    TERMINATE(isinf(matrix_get(M,j+dim,i)), "EM, i=%d,j=%d is inf", j+dim,i);
                    TERMINATE(isnan(matrix_get(M,j+dim,i)), "EM, i=%d,j=%d is nan", j+dim,i);
                }

                /* M_ME */
                {
                    sign_t sign;
                    edouble list_ij[] = { ln_bl1+cint.lnC_TM, ln_bl1+cint.lnD_TE };
                    sign_t signs_ij[] = { -sign_bl1*cint.signC_TM, -sign_bl1*cint.signD_TE};

                    edouble list_ji[] = { ln_bl2+cint.lnD_TM, ln_bl2+cint.lnC_TE };
                    sign_t signs_ji[] = { -MPOW(l1+l2+1)*sign_bl2*cint.signD_TM, -MPOW(l1+l2+1)*sign_bl2*cint.signC_TE };

                    matrix_set(M, i,dim+j, logadd_ms(list_ij, signs_ij, 2, &sign));
                    matrix_set(M_sign, i,dim+j, sign);

                    matrix_set(M, j,dim+i, logadd_ms(list_ji, signs_ji, 2, &sign));
                    matrix_set(M_sign, j,dim+i, sign);

                    TERMINATE(isinf(matrix_get(M,i,j+dim)), "ME, i=%d,j=%d is inf", i,j+dim);
                    TERMINATE(isnan(matrix_get(M,i,j+dim)), "ME, i=%d,j=%d is nan", i,j+dim);

                    TERMINATE(isinf(matrix_get(M,j,i+dim)), "ME, i=%d,j=%d is inf", j,i+dim);
                    TERMINATE(isnan(matrix_get(M,j,i+dim)), "ME, i=%d,j=%d is nan", j,i+dim);
                }
            }
        }
    }

    if(m == 0)
    {
        size_t i,j;
        matrix_edouble_t *EE = matrix_edouble_alloc(dim);
        matrix_edouble_t *MM = matrix_edouble_alloc(dim);

        matrix_sign_t *EE_sign = matrix_sign_alloc(dim);
        matrix_sign_t *MM_sign = matrix_sign_alloc(dim);

        for(i = 0; i < dim; i++)
            for(j = 0; j < dim; j++)
            {
                matrix_set(EE, i,j, matrix_get(M, i,j));
                matrix_set(EE_sign, i,j, matrix_get(M_sign, i,j));

                matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
                matrix_set(MM_sign, i,j, matrix_get(M_sign, i+dim,j+dim));
            }

        matrix_edouble_free(M);
        matrix_sign_free(M_sign);

        matrix_edouble_log_balance(MM);
        matrix_edouble_log_balance(EE);

        #ifdef USE_LAPACK
            logdet = matrix_logdet_lapack(EE, EE_sign) + matrix_logdet_lapack(MM, MM_sign);
        #else
            matrix_edouble_exp(EE, EE_sign);
            matrix_edouble_exp(MM, MM_sign);

            logdet = matrix_edouble_logdet(EE)+matrix_edouble_logdet(MM);
        #endif

        matrix_sign_free(EE_sign);
        matrix_sign_free(MM_sign);

        matrix_edouble_free(EE);
        matrix_edouble_free(MM);
    }
    else
    {
        matrix_edouble_log_balance(M);

        #ifdef USE_LAPACK
            logdet = matrix_logdet_lapack(M, M_sign);
        #else
            matrix_edouble_exp(M,M_sign);
            logdet = matrix_edouble_logdet(M);
        #endif

        matrix_edouble_free(M);
        matrix_sign_free(M_sign);
    }

    return logdet;
}

/*@}*/
