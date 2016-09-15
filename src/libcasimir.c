/**
 * @file   libcasimir.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   March, 2016
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */


/* for usleep and pthread_tryjoin_np */
#define _GNU_SOURCE

#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "floattypes.h"
#include "integration_drude.h"
#include "integration_perf.h"
#include "libcasimir.h"
#include "matrix.h"
#include "sfunc.h"
#include "utils.h"

/**
* @name Information on compilation and Casimir objects
*/
/*@{*/

/** @brief Return string with information about the binary
 *
 * The string contains date and time of compilation, the compiler and kind of
 * arithmetics the binary uses.
 *
 * This function is thread-safe.
 *
 * @param [out] str buffer for string
 * @param [in]  len length of string
 * @retval success bytes written if successful
 * @retval fail <0 if size == 0
 */
int casimir_compile_info(char *str, size_t size)
{
    if(size == 0)
        return -1;

    /* snprintf() writes at most size bytes (including the terminating null
     * byte ('\0')) to str. */
    return snprintf(str, size,
             "Compiled on %s at %s with %s",
              __DATE__, __TIME__, COMPILER
            );
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
    char buf[128];
    time_t timestamp = self->birthtime;
    struct tm ts = *localtime(&timestamp);
    strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", &ts);

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

    fprintf(stream, "%slmax            = %d\n", prefix, self->lmax);
    fprintf(stream, "%scores           = %d\n", prefix, self->cores);
    fprintf(stream, "%sprecision       = %g\n", prefix, self->precision);
    fprintf(stream, "%strace_threshold = %g\n", prefix, self->trace_threshold);
    fprintf(stream, "%sdetalg          = %s\n", prefix, self->detalg);
    fprintf(stream, "%sbalance         = %s\n", prefix, self->balance      ? "true" : "false");
    fprintf(stream, "%spivot           = %s\n", prefix, self->pivot        ? "true" : "false");
    fprintf(stream, "%sprecondition    = %s\n", prefix, self->precondition ? "true" : "false");
    fprintf(stream, "%scheck_elems     = %s\n", prefix, self->check_elems  ? "true" : "false");
    fprintf(stream, "%sbirthtime       = %s (%.1f)\n", prefix, buf, self->birthtime);
}


/**
 * @brief Thread safe wrapper to vfprintf
 *
 * This function is a wrapper to vfprintf and will print to stream, but it uses
 * locks to make this function thread-safe.
 *
 * @param [in] self Casimir object
 * @param [in] stream stream to print
 * @param [in] format format string
 * @param [in] args arguments for format string
 * @retval chars number of characters printed
 */
int casimir_vfprintf(casimir_t *self, FILE *stream, const char *format, va_list args)
{
    pthread_mutex_lock(&self->mutex);
    int ret = vfprintf(stream, format, args);
    pthread_mutex_unlock(&self->mutex);

    return ret;
}

/**
 * @brief Print debugging information
 *
 * Print debugging information to stderr if self->debug is set to true.
 *
 * A general note on debugging and error handling:
 *
 * 1) When an fatal error occures, for example allocating memory failed, the
 * program should terminate using the macro TERMINATE from utils.h.
 *
 * 2) If there was potentially a problem the program should use the macro WARN.
 * This macro will print a warning to stderr, but it will not terminate the
 * program. For example if there might be numerical instabilities or
 * convergence problems, but it is not 100% clear, the program should use WARN.
 *
 * 3) For debugging information / profiling information use casimir_debug.
 *
 * @param [in] self Casimir object
 * @param [in] format format string
 * @param [in] ... variables for for format string
 * @retval chars number of characters printed (see \ref casimir_vprintf)
 */
int casimir_debug(casimir_t *self, const char *format, ...)
{
    if(!self->debug)
        return 0;

    va_list args;
    va_start(args, format);
    int ret = casimir_vfprintf(self, stderr, format, args);
    va_end(args);

    return ret;
}

/**
 * @brief Print to stdout if flag verbose is set
 *
 * @param [in] self Casimir object
 * @param [in] format format string
 * @param [in] ... variables for for format string
 * @retval chars number of characters printed (see \ref casimir_vprintf)
 */
int casimir_verbose(casimir_t *self, const char *format, ...)
{
    if(!self->verbose)
        return 0;

    va_list args;
    va_start(args, format);
    int ret = casimir_vfprintf(self, stdout, format, args);
    va_end(args);

    return ret;
}

/*@}*/


/**
* @name various functions
*/
/*@{*/

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
float80 casimir_lnLambda(int l1, int l2, int m, sign_t *sign)
{
    if(sign != NULL)
        *sign = -1;
    return (log80(2*l1+1)+log80(2*l2+1)-log80(l1)-log80(l1+1)-log80(l2)-log80(l2+1)+lnfac80(l1-m)+lnfac80(l2-m)-lnfac80(l1+m)-lnfac80(l2+m))/2.0L;
}


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
float80 casimir_lnXi(int l1, int l2, int m, sign_t *sign)
{
    if(sign != NULL)
        *sign = MPOW(l2);
    return (log80(2*l1+1)+log80(2*l2+1)-lnfac80(l1-m)-lnfac80(l2-m)-lnfac80(l1+m)-lnfac80(l2+m)-log80(l1)-log80(l1+1)-log80(l2)-log80(l2+1))/2.0L \
           +lnfac80(2*l1)+lnfac80(2*l2)+lnfac80(l1+l2)-LOG4*(2*l1+l2+1)-lnfac80(l1-1)-lnfac80(l2-1);
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
    return 1+pow_2(omegap)/(xi*(xi+gamma_));
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
    return log1p(pow_2(omegap)/(xi*(xi+gamma_)));
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
void casimir_rp(casimir_t *self, float80 nT, float80 k, float80 *r_TE, float80 *r_TM)
{
    if(isinf(self->omegap_plane))
    {
        /* perfect reflectors */
        *r_TM = 1.0;
        *r_TE = -1.0;
    }
    else
    {
        /* Drude metals */
        const float80 epsilon = casimir_epsilon(nT, self->omegap_plane, self->gamma_plane);
        const float80 beta = sqrt80(1 + pow_2(nT)/(pow_2(nT)+pow_2(k)) * (epsilon-1));

        *r_TE = (1-beta)/(1+beta);
        *r_TM = (epsilon-beta)/(epsilon+beta);
    }
}

/*@}*/

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
* @name initialization and setting parameters
*/
/*@{*/

/**
 * @brief Create a new Casimir object for perfect reflectors
 *
 * This function will initialize a Casimir object with sphere and plane perfect
 * reflectors.
 *
 * By default, the value of lmax is chosen by:
 * lmax = ceil(max(CASIMIR_MINIMUM_LMAX, CASIMIR_FACTOR_LMAX/LbyR))
 *
 * Restrictions: \f$T > 0\f$, \f$0 < R/\mathcal{L} < 1\f$
 *
 * This function is not thread-safe.
 *
 * @param [out] self Casimir object
 * @param [in]  LbyR \f$\frac{L}{R}\f$
 * @param [in]  T temperature in units of \f$2\pi k_B \mathcal{L}/(\hbar c)\f$
 * @retval 0 if successful
 * @retval -1 if wrong value for LbyR
 * @retval -2 if wrong value for T
 */
int casimir_init(casimir_t *self, double LbyR, double T)
{
    TERMINATE(LDBL_MANT_DIG < 64, "No support for 80-bit extended precision, long double has %d bits.", LDBL_MANT_DIG);

    if(LbyR <= 0)
        return -1;
    if(T <= 0)
        return -2;

    self->lmax = ceil(MAX(CASIMIR_MINIMUM_LMAX, CASIMIR_FACTOR_LMAX/LbyR));

    self->T               = T;
    self->RbyScriptL      = 1./(1.+LbyR);
    self->LbyR            = LbyR;
    self->precision       = CASIMIR_PRECISION;
    self->trace_threshold = CASIMIR_TRACE_THRESHOLD;
    self->cores           = 1;
    self->threads         = xmalloc(self->cores*sizeof(pthread_t));

    /* initialize mie cache */
    casimir_mie_cache_init(self);

    /* perfect reflectors */
    self->omegap_sphere = INFINITY;
    self->gamma_sphere  = 0;
    self->omegap_plane  = INFINITY;
    self->gamma_plane   = 0;

    /* set verbose flag */
    self->verbose = false;

    /* initialize mutex for printf */
    pthread_mutex_init(&self->mutex, NULL);

    /**
     * parameters that users usually don't change
     */

    /* set debug flag */
    self->check_elems = true;

    /* set debug flag */
    self->debug = false;

    /* precondition matrix before QR decomposition */
    self->precondition = false;

    /* balance matrix before QR decomposition */
    self->balance = true;

    /* pivot matrix before QR decomposition */
    self->pivot = false;

    /* use QR decomposition to calculate determinant */
    memset(self->detalg, 0, sizeof(self->detalg));
    strcpy(self->detalg, CASIMIR_DETALG);

    self->birthtime = now();

    return 0;
}


/**
 * @brief Enable/disable debugging information
 *
 * This will enable/disable debugging information wether the flag debug is true
 * or false.
 *
 * @param [in] self Casimir object
 * @param [in] debug flag, true to enable, false to disable debugging
 */
void casimir_set_debug(casimir_t *self, bool debug)
{
    self->debug = debug;
}

bool casimir_get_debug(casimir_t *self)
{
    return self->debug;
}

void casimir_set_verbose(casimir_t *self, bool verbose)
{
    self->verbose = verbose;
}

bool casimir_get_verbose(casimir_t *self)
{
    return self->verbose;
}

/**
 * @brief Get birthtime
 *
 * The birthtime is the timestamp when the object was initialized.
 *
 * @param [in] self Casimir object
 * @retval timestamp
 */
double casimir_get_birthtime(casimir_t *self)
{
    return self->birthtime;
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
 * See \ref casimir_set_cores.
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
 * @brief Set algorithm to calculate deterimant
 *
 * The algorithm is given by detalg. Make sure that detalg contains a valid
 * algorithm, otherwise the computation will print a warning on runtime and
 * default to QR_FLOAT80.
 *
 * detalg may be: LU_FLOAT80, QR_FLOAT80, QR_LOG80 and (if supported)
 * QR_FLOAT128, QR_FLOATDD
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] detalg algorithm to compute deterimant
 * @retval 1
 */
int casimir_set_detalg(casimir_t *self, const char *detalg)
{
    strncpy(self->detalg, detalg, sizeof(self->detalg)/sizeof(char)-1);
    return 0;
}

/**
 * @brief Get algorithm to calculate deterimant
 *
 * The string is stored in detalg.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [out] detalg buffer
 * @retval 1
 */
int casimir_get_detalg(casimir_t *self, char detalg[128])
{
    strcpy(detalg, self->detalg);
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
 * @brief Get maximum value of \f$\ell\f$
 *
 * See \ref casimir_set_lmax.
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
 * @brief Get precision
 *
 * See \ref casimir_set_precision
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
 * @brief Set threshold for trace
 *
 * The threshold determines when Tr M is used as an approximation for
 * log(det(1-M)). If trace < threshold, then the value of the trace will be
 * used. Otherwise the determinant is caclulated.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] threshold threshold
 * @retval 1 if successful
 * @retval 0 if threshold < 0
 */
int casimir_set_trace_threshold(casimir_t *self, double threshold)
{
    if(threshold < 0)
        return 0;

    self->trace_threshold = threshold;
    return 1;
}

/**
 * @brief Get threshold for trace
 *
 * This function is not thread-safe.
 *
 * @param [in] self Casimir object
 * @retval threshold
 */
double casimir_get_trace_threshold(casimir_t *self)
{
    return self->trace_threshold;
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
    pthread_mutex_destroy(&self->mutex);

    casimir_mie_cache_free(self);

    xfree(self->threads);
    self->threads = NULL;
}

/*@}*/


/**
 * @name Mie coefficients
 */
/*@{*/

/**
 * @brief Return logarithm of prefactors \f$a_{\ell,0}^\mathrm{perf}\f$, \f$b_{\ell,0}^\mathrm{perf}\f$ and their signs
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
 * @param [in] l \f$\ell\f$
 * @param [out] a0 coefficient \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_a0 sign of \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] b0 coefficient \f$b_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_b0 sign of \f$b_{\ell,0}^\mathrm{perf}\f$
 */
void casimir_lnab0(int l, float80 *a0, sign_t *sign_a0, float80 *b0, sign_t *sign_b0)
{
    *sign_a0 = MPOW(l);
    *sign_b0 = MPOW(l+1);
    *b0 = LOGPI-lngamma80(l+0.5)-lngamma80(l+1.5);
    *a0 = *b0+log1p80(1.0L/l);
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
 * @param [in] l \f$\ell\f$
 * @param [in] n Matsubara term, \f$xi = nT\f$
 * @param [out] sign sign of \f$a_\ell\f$
 * @retval logarithm of Mie coefficient \f$a_\ell\f$
 */
double casimir_lna_perf(casimir_t *self, const int l, const int n, sign_t *sign)
{
    float80 lnKlp,lnKlm,lnIlm,lnIlp;
    const float80 chi = n*self->T*self->RbyScriptL;

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother?
     */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp);

    /* We want to calculate
     * a(chi) = (-1)^(l+1)*pi/2 * ( l*Ip-chi*Im )/( l*Kp+chi*Km )
     *        = (-1)^(l+1)*pi/2*Ip/Kp * ( l-chi*Im/Ip )/( l+chi*Km/Kp )
     *          \--------/ \--------/   \-------------/ \-------------/
     *             sign     prefactor      numerator      denominator
     */

    const float80 prefactor = LOGPI-LOG2+lnIlp-lnKlp;

    /* numerator */
    const float80 numerator = logadd_s(log80(l), +1, log80(chi)+lnIlm-lnIlp, -1, sign);
    /* denominator */
    const float80 denominator = logadd(log80(l), log80(chi)+lnKlm-lnKlp);

    *sign *= MPOW(l+1);
    return prefactor+numerator-denominator;
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
 * @param [in] l \f$\ell\f$
 * @param [in] n Matsubara term, \f$\xi = nT\f$
 * @param [out] sign sign of \f$b_\ell\f$
 * @retval logarithm of Mie coefficient \f$b_\ell\f$
 */
double casimir_lnb_perf(casimir_t *self, const int l, const int n, sign_t *sign)
{
    const float80 chi = n*self->T*self->RbyScriptL;
    float80 lnInu, lnKnu;

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
 * sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * For Drude metals we calculate the Mie coefficients al(iξ) und bl(iξ) using
 * the expressions taken from [1]. Ref. [1] is the erratum to [2]. Please note
 * that the equations (3.30) and (3.31) in [3] are wrong.
 *
 * Note: If sla =~ slb or slc =~ sld, there is a loss of significance when
 * calculating sla-slb or slc-sld.
 *
 * This function is thread safe - as long you don't change temperature, aspect
 * ratio and dielectric properties of sphere.
 *
 * References:
 * [1] Erratum: Thermal Casimir effect for Drude metals in the plane-sphere
 * geometry, Canaguier-Durand, Neto, Lambrecht, Reynaud (2010)
 * http://journals.aps.org/pra/abstract/10.1103/PhysRevA.83.039905
 * [2] Thermal Casimir effect for Drude metals in the plane-sphere geometry,
 * Canaguier-Durand, Neto, Lambrecht, Reynaud (2010),
 * http://journals.aps.org/pra/abstract/10.1103/PhysRevA.82.012511
 * [3] Negative Casimir entropies in the plane-sphere geometry, Hartmann, 2014
 *
 * @param [in,out] self Casimir object
 * @param [in] n_mat Matsubara term, \f$\xi = nT\f$
 * @param [in] l \f$\ell\f$
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 * @param [out] sign_a sign of Mie coefficient \f$a_\ell\f$
 * @param [out] sign_b sign of Mie coefficient \f$b_\ell\f$
 */
void casimir_lnab(casimir_t *self, const int n_mat, const int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
{
    if(isinf(self->omegap_sphere))
    {
        /* Mie coefficients for perfect reflectors */
        *lna = casimir_lna_perf(self, l, n_mat, sign_a);
        *lnb = casimir_lnb_perf(self, l, n_mat, sign_b);
        return;
    }

    /* Mie coefficients for Drude metals */

    /* ξ = nT */
    const float80 xi = n_mat*self->T;

    /* χ = ξ*R/(R+L) = ξ/(1+L/R) */
    const float80 chi    = xi/(1+self->LbyR);
    const float80 ln_chi = log80(xi)-log1p80(self->LbyR);
    const float80 ln_l   = log80(l);

    /**
     * Note: n is the refraction index, n_mat the Matsubara index
     * n    = sqrt(ε(ξ,ω_p,γ))
     * ln_n = ln(sqrt(ε(ξ,ω_p,γ)))
     */
    const float80 ln_n = casimir_lnepsilon(xi, self->omegap_sphere, self->gamma_sphere)/2;
    const float80 n    = exp80(ln_n);

    float80 lnIlp, lnKlp, lnIlm, lnKlm, lnIlp_nchi, lnKlp_nchi, lnIlm_nchi, lnKlm_nchi;

    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp); /* I_{l+0.5}(χ), K_{l+0.5}(χ) */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm); /* K_{l-0.5}(χ), K_{l-0.5}(χ) */

    bessel_lnInuKnu(l,   n*chi, &lnIlp_nchi, &lnKlp_nchi); /* I_{l+0.5}(nχ), K_{l+0.5}(nχ) */
    bessel_lnInuKnu(l-1, n*chi, &lnIlm_nchi, &lnKlm_nchi); /* K_{l-0.5}(nχ), K_{l-0.5}(nχ) */

    sign_t sign_sla, sign_slb, sign_slc, sign_sld;
    const float80 ln_sla = lnIlp_nchi + logadd_s(ln_l+lnIlp,      +1,      ln_chi+lnIlm,      -1, &sign_sla);
    const float80 ln_slb = lnIlp      + logadd_s(ln_l+lnIlp_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_slb);
    const float80 ln_slc = lnIlp_nchi + logadd_s(ln_l+lnKlp,      +1,      ln_chi+lnKlm,      +1, &sign_slc);
    const float80 ln_sld = lnKlp      + logadd_s(ln_l+lnIlp_nchi, +1, ln_n+ln_chi+lnIlm_nchi, -1, &sign_sld);

    /**
     */
    /*
    printf("n =%.18Lg\n",     exp80(ln_n));
    printf("n2=%.18Lg\n",     exp80(2*ln_n));
    printf("lnIl = %.18Lg\n", exp80(lnIl));
    printf("chi=%.18Lg\n",    chi);

    printf("sla=%.18Lg\n", sign_sla*exp80(ln_sla));
    printf("slb=%.18Lg\n", sign_slb*exp80(ln_slb));
    printf("slc=%.18Lg\n", sign_slc*exp80(ln_slc));
    printf("sld=%.18Lg\n", sign_sld*exp80(ln_sld));
    */

    sign_t sign_a_num, sign_a_denom, sign_b_num, sign_b_denom;
    *lna = LOGPI - LOG2 + logadd_s(2*ln_n+ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_a_num) - logadd_s(2*ln_n+ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_a_denom);
    *lnb = LOGPI - LOG2 + logadd_s(       ln_sla, +sign_sla, ln_slb, -sign_slb, &sign_b_num) - logadd_s(       ln_slc, +sign_slc, ln_sld, -sign_sld, &sign_b_denom);

    *sign_a = MPOW(l+1)*sign_a_num*sign_a_denom;
    *sign_b = MPOW(l+1)*sign_b_num*sign_b_denom;
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
 * @param [in,out] self Casimir object
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
 * You usually don't want to use this function yourself.
 *
 * This function is not thread safe.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term
 */
void casimir_mie_cache_alloc(casimir_t *self, int n)
{
    const int nmax = self->mie_cache->nmax;
    const int lmax = self->mie_cache->lmax;
    casimir_mie_cache_t *cache = self->mie_cache;

    if(n > nmax)
    {
        cache->entries = xrealloc(cache->entries, (n+1)*sizeof(casimir_mie_cache_entry_t *));

        for(int l = nmax+1; l <= n; l++)
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
        for(int l = 1; l <= lmax; l++)
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
 * @param [in, out] self Casimir object
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
 * @param [in, out] self Casimir object
 * @param [in] l \f$\ell\f$
 * @param [in] n Matsubara term
 * @param [out] ln_a logarithm of \f$a_\ell\f$
 * @param [out] sign_a sign of \f$a_\ell\f$
 * @param [out] ln_b logarithm of \f$b_\ell\f$
 * @param [out] sign_b sign of \f$b_\ell\f$
 */
void casimir_mie_cache_get(casimir_t *self, int l, int n, double *ln_a, sign_t *sign_a, double *ln_b, sign_t *sign_b)
{
    /* this mutex is important to prevent memory corruption */
    pthread_mutex_lock(&self->mie_cache->mutex);

    int nmax = self->mie_cache->nmax;
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
    casimir_mie_cache_entry_t *entry = self->mie_cache->entries[n];

    /* at this point it is finally is safe to release the mutex */
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
 * @param [in, out] self Casimir object
 */
void casimir_mie_cache_free(casimir_t *self)
{
    casimir_mie_cache_t *cache = self->mie_cache;
    casimir_mie_cache_entry_t **entries = cache->entries;

    pthread_mutex_destroy(&cache->mutex);

    /* free
     * 1) the lists of al, bl, sign_al, sign_bl for every entry
     * 2) every entry (casimir_mie_cache_entry_t)
     * 3) the list of entries
     * 4) the mie cache object (casimir_mie_cache_t)
     */
    for(int n = 0; n <= cache->nmax; n++)
    {
        if(entries[n] != NULL)
        {
            xfree(entries[n]->ln_al);
            xfree(entries[n]->sign_al);
            xfree(entries[n]->ln_bl);
            xfree(entries[n]->sign_bl);
            xfree(entries[n]);

            entries[n] = NULL;
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
static double _sum(double values[], int len)
{
    double sum = 0;

    for(int i = len-1; i > 0; i--)
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
    int joined = 0, running = 0;
    casimir_thread_t *r;
    pthread_t **threads = self->threads;

    for(int i = 0; i < self->cores; i++)
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
 * ratio, dielectric properties of sphere and plane, lmax and integration.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term, \f$\xi=nT\f$
 * @param [out] mmax maximum number of \f$m\f$
 * @retval F Casimir free energy for given \f$n\f$
 */
double casimir_F_n(casimir_t *self, const int n, int *mmax)
{
    int m;
    const double precision = self->precision;
    const int lmax = self->lmax;
    double sum_n = 0;
    double values[lmax+1];

    /* initialize to 0 */
    for(m = 0; m <= lmax; m++)
        values[m] = 0;

    for(m = 0; m <= self->lmax; m++)
    {
        values[m] = casimir_logdetD(self,n,m);
        casimir_verbose(self, "# n=%d, m=%d, logdetD=%.15g\n", n, m, values[m]);

        /* The stop criterion is as follows: If
         * logdetD(n,m)/(\sum_i^m logdetD(n,i)) < precision
         * we can skip the summation.
         */
        sum_n = _sum(values, lmax+1);
        if(values[0] != 0 && fabs(values[m]/sum_n) < precision)
            break;
    }

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
    int n = 0;
    double sum_n = 0;
    const double precision = self->precision;
    double *values = NULL;
    int len = 0;
    int ncalc = 0;
    const int cores = self->cores;
    const int delta = MAX(1024, cores);
    pthread_t **threads = self->threads;

    for(int i = 0; i < cores; i++)
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

            for(int i = len; i < len+delta; i++)
                values[i] = 0;

            len += delta;
        }

        if(cores > 1)
        {
            for(int i = 0; i < cores; i++)
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
    const float80 lnRbyScriptL = log80(self->RbyScriptL);
    matrix_float80 *EE = NULL, *MM = NULL;
    matrix_sign_t *EE_sign = NULL, *MM_sign = NULL;

    const int min = MAX(m,1);
    const int max = self->lmax;
    const int dim = (max-min+1);

    if(logdet_EE != NULL)
    {
        EE = matrix_float80_alloc(dim);
        EE_sign = matrix_sign_alloc(dim);
    }
    if(logdet_MM != NULL)
    {
        MM = matrix_float80_alloc(dim);
        MM_sign = matrix_sign_alloc(dim);
    }

    /* calculate the logarithm of the matrix elements of -M. The function
     * matrix_logdetIdpM then calculates log(det(1-M)) = log(det(D)) */
    for(int l1 = min; l1 <= max; l1++)
    {
        sign_t sign_a0, sign_b0;
        float80 lna0, lnb0;
        casimir_lnab0(l1, &lna0, &sign_a0, &lnb0, &sign_b0);

        for(int l2 = min; l2 <= max; l2++)
        {
            /* i: row of matrix, j: column of matrix */
            const int i = l1-min, j = l2-min;
            sign_t sign_xi;
            const float80 lnXiRL = casimir_lnXi(l1,l2,m,&sign_xi)+(2*l1+1)*lnRbyScriptL;

            if(EE != NULL)
            {
                matrix_set(EE, i,j, lna0+lnXiRL);
                matrix_set(EE_sign, i,j, -sign_xi*sign_a0);
            }
            if(MM != NULL)
            {
                matrix_set(MM, i,j, lnb0+lnXiRL);
                matrix_set(MM_sign, i,j, sign_xi*sign_b0);
            }
        }
    }

    /* calculate logdet and free space */
    if(EE != NULL)
    {
        *logdet_EE = matrix_logdetIdpM(self, EE, EE_sign);

        matrix_sign_free(EE_sign);
        matrix_float80_free(EE);
    }
    if(MM != NULL)
    {
        *logdet_MM = matrix_logdetIdpM(self, MM, MM_sign);

        matrix_sign_free(MM_sign);
        matrix_float80_free(MM);
    }
}


/**
 * @brief Calculate \f$\mathrm{tr} \mathcal{D}^{(m)}(\xi=nT)\f$
 *
 * This function calculates the trace of the scattering operator D for the
 * Matsubara term \f$n\f$.
 *
 * This function is thread-safe - as long you don't change lmax, temperature,
 * aspect ratio, dielectric properties of sphere and plane, and integration.
 *
 * @param [in,out] self Casimir object
 * @param [in] n Matsubara term
 * @param [in] m
 * @param [in,out] integration_obj may be NULL
 * @retval trM \f$\mathrm{tr} \mathcal{D}^{(m)}(\xi=nT)\f$
 */
double casimir_trM(casimir_t *self, int n, int m, void *obj)
{
    const int min = MAX(m,1);
    const int max = self->lmax;
    integration_perf_t  *int_perf  = obj;
    integration_drude_t *int_drude = obj;
    float80 trM = 0;

    for(int l = min; l <= max; l++)
    {
        float80 sum;
        casimir_integrals_t cint;
        double ln_al, ln_bl;
        sign_t sign, sign_al, sign_bl;

        casimir_mie_cache_get(self, l, n, &ln_al, &sign_al, &ln_bl, &sign_bl);

        if(isinf(self->omegap_plane))
            casimir_integrate_perf(int_perf, l, l, &cint);
        else
            casimir_integrate_drude(int_drude, l, l, &cint);

        /* EE */
        /* A_TE + B_TM */
        sum = logadd_s(cint.lnA_TE, cint.signA_TE, cint.lnB_TM, cint.signB_TM, &sign);
        trM += sign*sign_al*exp80(ln_al+sum);

        /* MM */
        /* A_TM + B_TE */
        sum = logadd_s(cint.lnA_TM, cint.signA_TM, cint.lnB_TE, cint.signB_TE, &sign);
        trM += sign*sign_bl*exp80(ln_bl+sum);
    }

    return trM;
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
 * @retval logdetD \f$\log \det \mathcal{D}^{(m)}(\xi=nT)\f$
 */
double casimir_logdetD(casimir_t *self, int n, int m)
{
    const double start = now();
    const double nT = n*self->T;
    double logdet = 0;
    integration_perf_t int_perf;
    
    integration_drude_t int_drude = {
        .plm_cache = NULL,
        .m         = m,
        .nT        = n * self->T
    };

    const int min = MAX(m,1);
    const int max = self->lmax;
    const int dim = (max-min+1);

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

    double trace;
    if(isinf(self->omegap_plane))
    {
        /* perfect reflector */
        casimir_integrate_perf_init(&int_perf, nT, m, self->lmax);
        trace = casimir_trM(self, n, m, &int_perf);
    }
    else
    {
        /* drude mirror */
        casimir_integrate_drude_init(self, &int_drude, nT, m, self->lmax);
        trace = casimir_trM(self, n, m, &int_drude);
    }

    /* If |tr M| is smaller than trace_threshold, we use the approximation
     *      log det D = log det(Id-M) ≈ -tr M
     * instead of calculating all matrix elements of M and computing the
     * determinant. So, if the trace of the round trip matrix M is small, we
     * don't have to compute O(lmax²), but only O(lmax) matrix elements. This
     * will drastically speed up the calculation!
     */
    if(fabs(trace) < self->trace_threshold)
    {
        /* free integration object */
        if(isinf(self->omegap_plane))
            casimir_integrate_perf_free(&int_perf);
        else
            casimir_integrate_drude_free(&int_drude);

        casimir_debug(self, "# calculated %dx%d matrix elements (trace approximation): %gs\n", 2*dim, 2*dim, now()-start);

        return -trace;
    }

    matrix_float80 *M     = matrix_float80_alloc(2*dim);
    matrix_sign_t *M_sign = matrix_sign_alloc   (2*dim);

    /* set matrix elements to NAN */
    if(self->check_elems)
    {
        for(int i = 0; i < pow_2(M->dim); i++)
            M->M[i] = NAN;
    }

    /* M_EE, -M_EM
       M_ME,  M_MM */
    /*
     * Theese two for-loops are a little bit difficult to understand.
     * They do the same as
     *     for(int l1 = min; l1 <= max; ++l1)
     *         for(int l2 = min; l2 <= l1; ++l2).
     * But they change the order of the evaluation.
     * This doesn't change the result, since the order of evaluation doesn't matter. But it is faster
     * for Drude, because we can reuse most of the legendre polynomials (see plm_cache.c).
     * This will successively evaluate the matrix-elements with the same "l1 + l2".
     */
    for(int l1_plus_l2 = 2 * min; l1_plus_l2 <= 2 * max; ++l1_plus_l2)
    {
        for(int l2 = min; l2 <= l1_plus_l2 / 2; ++l2)
        {
            const int l1 = l1_plus_l2 - l2;
            if(l1 > max)
                continue;

            const int i = l1-min;
            const int j = l2-min;

            double ln_al1, ln_bl1, ln_al2, ln_bl2;
            sign_t sign_al1, sign_bl1, sign_al2, sign_bl2;

            casimir_mie_cache_get(self, l1, n, &ln_al1, &sign_al1, &ln_bl1, &sign_bl1);
            casimir_mie_cache_get(self, l2, n, &ln_al2, &sign_al2, &ln_bl2, &sign_bl2);

            casimir_integrals_t cint;
            if(isinf(self->omegap_plane))
                casimir_integrate_perf(&int_perf, l1, l2, &cint);
            else
                casimir_integrate_drude(&int_drude, l1, l2, &cint);

            /* EE */
            {
                sign_t sign;
                /* A_TE + B_TM */
                float80 sum = logadd_s(cint.lnA_TE, cint.signA_TE, cint.lnB_TM, cint.signB_TM, &sign);

                matrix_set(M, i,j, ln_al1+sum);
                matrix_set(M, j,i, ln_al2+sum);

                matrix_set(M_sign, i,j, -            sign_al1*sign);
                matrix_set(M_sign, j,i, -MPOW(l1+l2)*sign_al2*sign);
            }

            /* MM */
            {
                sign_t sign;
                /* A_TM + B_TE */
                float80 sum = logadd_s(cint.lnA_TM, cint.signA_TM, cint.lnB_TE, cint.signB_TE, &sign);

                matrix_set(M, i+dim,j+dim, ln_bl1+sum);
                matrix_set(M, j+dim,i+dim, ln_bl2+sum);

                matrix_set(M_sign, i+dim,j+dim, -            sign_bl1*sign);
                matrix_set(M_sign, j+dim,i+dim, -MPOW(l1+l2)*sign_bl2*sign);
            }


            if(m != 0)
            {
                sign_t sign1, sign2;

                /* C_TE + D_TM */
                float80 sum1 = logadd_s(cint.lnC_TE, cint.signC_TE, cint.lnD_TM, cint.signD_TM, &sign1);

                /* C_TM + D_TE */
                float80 sum2 = logadd_s(cint.lnC_TM, cint.signC_TM, cint.lnD_TE, cint.signD_TE, &sign2);

                /* M_EM */
                {
                    matrix_set(M, i,dim+j, ln_al1+sum1);
                    matrix_set(M, j,dim+i, ln_al2+sum2);

                    matrix_set(M_sign, i,dim+j, -              sign_al1*sign1);
                    matrix_set(M_sign, j,dim+i, -MPOW(l1+l2+1)*sign_al2*sign2);
                }

                /* M_ME */
                {
                    matrix_set(M, dim+i,j, ln_bl1+sum2);
                    matrix_set(M, dim+j,i, ln_bl2+sum1);

                    matrix_set(M_sign, dim+i,j, -              sign_bl1*sign2);
                    matrix_set(M_sign, dim+j,i, -MPOW(l1+l2+1)*sign_bl2*sign1);
                }
            }
        }
    }

    if(isinf(self->omegap_plane))
        casimir_integrate_perf_free(&int_perf);
    else
        casimir_integrate_drude_free(&int_drude);

    casimir_debug(self, "# calculating %dx%d matrix elements (trace approximation): %gs\n", 2*dim, 2*dim, now()-start);

    /* check if matrix elements are finite */
    if(self->check_elems)
    {
        for(int l1 = min; l1 <= max; l1++)
            for(int l2 = min; l2 <= max; l2++)
            {
                const int i = l1-min;
                const int j = l2-min;
                const float80 elem_EE = matrix_get(M, i,j);
                const float80 elem_MM = matrix_get(M, i+dim,j+dim);

                TERMINATE(!isfinite(elem_EE), "matrix element not finite: P=EE, l1=%d, l2=%d, elem=%Lg", l1,l2, elem_EE);
                TERMINATE(!isfinite(elem_MM), "matrix element not finite: P=MM, l1=%d, l2=%d, elem=%Lg", l1,l2, elem_MM);

                if(m != 0)
                {
                    const float80 elem_EM = matrix_get(M, i+dim,j);
                    const float80 elem_ME = matrix_get(M, i,j+dim);

                    TERMINATE(!isfinite(elem_EM), "matrix element not finite: P=EM, l1=%d, l2=%d, elem=%Lg", l1,l2, elem_EM);
                    TERMINATE(!isfinite(elem_ME), "matrix element not finite: P=EM, l1=%d, l2=%d, elem=%Lg", l1,l2, elem_EM);
                }

            }
    }

#if 0
    /* Dump matrix */
    int ret;

    ret = matrix_float80_save(M, "M.out");
    WARN(ret, "Couldn't dump matrix M.out");

    ret = matrix_sign_save(M_sign, "M_sign.out");
    WARN(ret, "Couldn't dump matrix M_sign.out");
#endif

    /* We have calculated -M here. We now call matrix_logdetIdpM that will
     * calculate log(det(1-M)) = log(det(D)) */

    if(m == 0)
    {
        /* The memory footprint can be improved here */
        matrix_float80 *EE = matrix_float80_alloc(dim);
        matrix_float80 *MM = matrix_float80_alloc(dim);

        matrix_sign_t *EE_sign = matrix_sign_alloc(dim);
        matrix_sign_t *MM_sign = matrix_sign_alloc(dim);

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
            {
                matrix_set(EE, i,j, matrix_get(M, i,j));
                matrix_set(EE_sign, i,j, matrix_get(M_sign, i,j));

                matrix_set(MM, i,j, matrix_get(M, dim+i,dim+j));
                matrix_set(MM_sign, i,j, matrix_get(M_sign, i+dim,j+dim));
            }

        matrix_float80_free(M);
        matrix_sign_free(M_sign);

        logdet  = matrix_logdetIdpM(self, EE, EE_sign);
        logdet += matrix_logdetIdpM(self, MM, MM_sign);

        matrix_sign_free(EE_sign);
        matrix_sign_free(MM_sign);

        matrix_float80_free(EE);
        matrix_float80_free(MM);
    }
    else
    {
        logdet = matrix_logdetIdpM(self, M, M_sign);

        matrix_float80_free(M);
        matrix_sign_free(M_sign);
    }

    return logdet;
}

/*@}*/
