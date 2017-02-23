/**
 * @file   libcasimir.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   February, 2017
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */

#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <hodlr.h>

#include "integration.h"
#include "libcasimir.h"
#include "matrix.h"
#include "sfunc.h"
#include "utils.h"

/**
* @name Information on compilation and Casimir objects
*/
/*@{*/

/** @brief Return string with information about compilation
 *
 * The string contains date and time of compilation, and the name of the
 * compiler.
 *
 * @param [out] str buffer for string
 * @param [in]  len length of string
 * @retval success bytes written if successful
 */
int casimir_compile_info(char *str, size_t size)
{
    /* snprintf() writes at most size bytes (including the terminating null
     * byte \0) to str. */
    return snprintf(str, size, "Compiled on %s at %s with %s", __DATE__, __TIME__, COMPILER);
}


/** @brief Print object information to stream
 *
 * Print information about the object self to stream.
 *
 * @param self Casimir object
 * @param stream where to print the string
 * @param prefix if prefix != NULL: start every line with the string contained
 * in prefix
 */
void casimir_info(casimir_t *self, FILE *stream, const char *prefix)
{
    const char *detalg_str;
    if(prefix == NULL)
        prefix = "";

    switch(self->detalg)
    {
        case DETALG_HODLR: detalg_str = "HODLR"; break;
        case DETALG_LU:    detalg_str = "LU";    break;
        default:           detalg_str = "unknown";
    }

    fprintf(stream, "%sL/R       = %.8g\n", prefix, self->LbyR);
    fprintf(stream, "%sldim      = %d\n", prefix, self->ldim);
    fprintf(stream, "%stolerance = %g\n", prefix, self->tolerance);
    fprintf(stream, "%sthreshold = %g\n", prefix, self->threshold);
    fprintf(stream, "%sdetalg    = %s\n", prefix, detalg_str);
    fprintf(stream, "%sverbose   = %s\n", prefix, self->verbose ? "true" : "false");
    fprintf(stream, "%sdebug     = %s\n", prefix, self->debug ? "true" : "false");
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
 * 1) When a fatal error occures, e.g. allocating memory failed, the program
 * should terminate using the macro TERMINATE from utils.h.
 *
 * 2) If there was probably a problem the program should use the macro WARN.
 * This macro will print a warning to stderr, but it will not terminate the
 * program. For example, if there might be numerical instabilities or
 * convergence problems, the program should use WARN.
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
 * @brief Calculate logarithm \f$-\Lambda_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of \f$-\Lambda_{\ell_1 \ell_2}^{(m)}\f$ for
 * \f$\ell_1,\ell_2,m\f$. This prefactor is defined by (cf Eq. (5.19))
 * \f[
 *      \Lambda_{\ell_1,\ell_2}^{(m)} = -\frac{2 N_{\ell_1,m} N_{\ell_2,m}}{\sqrt{\ell_1 (\ell_1+1) \ell_2 (\ell_2+1)}}
 * \f]
 *
 * The values are computed using the lgamma function to avoid overflows.
 *
 * Symmetries: \f$\Lambda_{\ell_1,\ell_2}^{(m)} = \Lambda_{\ell_2,\ell_1}^{(m)}\f$
 *
 * @param [in]  l1 \f$\ell_1\f$, l1>0
 * @param [in]  l2 \f$\ell_2\f$, l2>0
 * @param [in]  m  \f$m\f$, MIN(l1,l2)>=m>=0
 * @retval lnLambda \f$\log{\Lambda_{\ell_1,\ell_2}^{(m)}}\f$
 */
double casimir_lnLambda(int l1, int l2, int m)
{
    return (logi(2*l1+1)+logi(2*l2+1)-logi(l1)-logi(l1+1)-logi(l2)-logi(l2+1)+lfac(l1-m)+lfac(l2-m)-lfac(l1+m)-lfac(l2+m))/2.0;
}

/**
 * @brief Evaluate dielectric function
 *
 * The dielectrict function set by casimir_set_epsilonm1 will be called.
 *
 * @param [in] self Casimir object
 * @param [in] xi
 * @retval epsilon-1, epsilon(xi)-1
 */
double casimir_epsilonm1(casimir_t *self, double xi)
{
    return self->epsilonm1(xi, self->userdata);
}

/**
 * @brief Dielectric function for perfect reflectors
 *
 * @param [in] xi ignored
 * @param [in] userdata ignored
 * @retval INFINITY epsilon(xi) = INFINITY
 */
double casimir_epsilonm1_perf(__attribute__((unused)) double xi, __attribute__((unused)) void *userdata)
{
    return INFINITY;
}

/**
 * @brief Dielectric function for Drude reflectors
 *
 * Dielectric function for Drude
 *      epsilon(ξ) = ωp²/(ξ*(ξ+γ))
 *
 * The parameters ωp and γ must be provided by userdata:
 *      omegap = userdata[0]
 *      gamma_ = userdata[1]
 *
 * @param [in] xi dielectric function
 * @param [in] userdata userdata
 * @retval epsilon epsilon(xi)
 */
double casimir_epsilonm1_drude(double xi, void *userdata)
{
    double *ptr = (double *)userdata;
    double omegap = ptr[0], gamma_ = ptr[1];

    return pow_2(omegap)/(xi*(xi+gamma_));
}

/**
 * @brief Calculate Fresnel coefficients \f$r_{TE}\f$ and \f$r_{TM}\f$ for arbitrary metals
 *
 * This function calculates the Fresnel coefficients for TE and TM mode
 *
 * @param [in]     self  Casimir object
 * @param [in]     nT    \f$\xi=nT\f$ imaginary frequency
 * @param [in]     k     xy projection of wavevector
 * @param [in,out] r_TE  Fresnel coefficient for TE mode
 * @param [in,out] r_TM  Fresnel coefficient for TM mode
 */
void casimir_rp(casimir_t *self, double nT, double k, double *r_TE, double *r_TM)
{
    const double epsilonm1 = casimir_epsilonm1(self, nT);

    if(isinf(epsilonm1))
    {
        /* perfect reflectors */
        *r_TM = 1.0;
        *r_TE = -1.0;

        /* get us out of here */
        return;
    }

    /* Arbitrary metals
     *
     * In scaled units
     *     β = sqrt( 1 + ξ²/(ξ²+k²)*(ε-1) ) = sqrt(1+x),
     * where
     *     x = ξ²/(ξ²+k²)*(ε-1).
     *
     * We calculate x. If x is small, β≈1 and a loss of significance
     * occures when calculating 1-β.
     *
     * For this reason we use the Taylor series
     *     sqrt(1+x) ≈ 1 + x/2 - x²/8 + x³/16 - 5*x^4/128 + ...
     * to avoid a loss of significance if x is small.
     *
     * Note: ξ=nT
     */
    const double x = pow_2(nT)/(pow_2(nT)+pow_2(k))*epsilonm1;

    if(fabs(x) < 1e-5)
    {
        /* β-1 = sqrt(1+x)-1 = x/2 - x²/8 + x³/16 - 5*x^4/128 + O(x^5) */
        const double betam1 = x/2 - pow_2(x)/8 + pow_3(x)/16 - 5*pow_4(x)/128;

        *r_TE = -betam1/(2+betam1);
        *r_TM = (epsilonm1-betam1)/(epsilonm1+2+betam1);
    }
    else
    {
        const double beta = sqrt(1+x);

        *r_TE = (1-beta)/(1+beta);
        *r_TM = (epsilonm1+1-beta)/(epsilonm1+1+beta);
    }
}

/*@}*/

/**
* @name initialization and setting parameters
*/
/*@{*/

/**
 * @brief Create a new Casimir object
 *
 * This function will initialize a Casimir object.  By default the dielectric
 * function corresponds to perfect reflectors, i.e. epsilon=inf.
 *
 * By default, the value of ldim is chosen by:
 * ldim = ceil(max(CASIMIR_MINIMUM_LDIM, CASIMIR_FACTOR_LDIM/LbyR))
 *
 * Restrictions: \f$\frac{L}{R} > 0\f$
 *
 * @param [out] self Casimir object
 * @param [in]  LbyR \f$\frac{L}{R}\f$, LbyR > 0
 * @retval object Casimir object if successful
 * @retval NULL   an error occured
 */
casimir_t *casimir_init(double LbyR)
{
    if(LbyR <= 0)
        return NULL;

    casimir_t *self = xmalloc(sizeof(casimir_t));

    /* geometry */
    self->LbyR = LbyR;
    self->y = -M_LOG2-log1p(LbyR);

    /* dimension of vector space */
    self->ldim = ceil(MAX(CASIMIR_MINIMUM_LDIM, CASIMIR_FACTOR_LDIM/LbyR));

    /* perfect reflectors */
    self->epsilonm1 = casimir_epsilonm1_perf;
    self->userdata  = NULL;

    /* set threshold */
    self->threshold = CASIMIR_THRESHOLD;

    /* integration error */
    self->tolerance = CASIMIR_TOLERANCE;

    /* set verbose flag */
    self->verbose = false;

    /* initialize mutex for printf */
    pthread_mutex_init(&self->mutex, NULL);

    /* set debug flag */
    self->debug = false;

    /* use LU decomposition by default */
    self->detalg = DETALG_HODLR;

    return self;
}


/**
 * @brief Free memory for Casimir object
 *
 * Free allocated memory for the Casimir object self.
 *
 * @param [in,out] self Casimir object
 */
void casimir_free(casimir_t *self)
{
    if(self != NULL)
    {
        pthread_mutex_destroy(&self->mutex);
        xfree(self);
    }
}

/**
 * @brief Enable/disable debugging information
 *
 * Enable/disable debugging information if flag debug is true/false.
 *
 * @param [in] self Casimir object
 * @param [in] debug flag, true to enable, false to disable debugging
 */
void casimir_set_debug(casimir_t *self, bool debug)
{
    self->debug = debug;
}

/**
 * @brief Get debugging flag
 *
 * @retval debug true/false
 */
bool casimir_get_debug(casimir_t *self)
{
    return self->debug;
}

/**
 * @brief Enable/disable verbose information
 *
 * Enable/disable verbose information if flag debug is true/false.
 *
 * @param [in] self Casimir object
 * @param [in] verbose flag, true to enable, false to disable
 */
void casimir_set_verbose(casimir_t *self, bool verbose)
{
    self->verbose = verbose;
}

/**
 * @brief Get verbose flag
 *
 * @retval verbose true/false
 */
bool casimir_get_verbose(casimir_t *self)
{
    return self->verbose;
}

/**
 * @brief Set tolerance for numerical integration
 *
 * Set relative tolerance for numerical integration.
 *
 * @param [in] self Casimir object
 * @param [in] tolerance relative tolerance
 * @retval 0 if an error occured
 * @retval 1 on success
 */
int casimir_set_tolerance(casimir_t *self, double tolerance)
{
    if(tolerance <= 0)
        return 0;

    self->tolerance = tolerance;
    return 1;
}

/**
 * @brief Get relative tolerance for numerical integration
 *
 * @retval tolerance numerical tolerance
 */
double casimir_get_tolerance(casimir_t *self)
{
    return self->tolerance;
}

/**
 * @brief Set material parameters for generic metals
 *
 * Set dielectric function of material.
 *
 * If epsilonm1 is NULL, assume perfect reflectors.
 *
 * @param [in,out] self Casimir object
 * @param [in] epsilonm1  callback to the function that calculates epsilon(i*xi)-1
 * @param [in] userdata   arbitrary pointer to data that is passwd to epsilonm1 whenever the function is called
 */
void casimir_set_epsilonm1(casimir_t *self, double (*epsilonm1)(double xi, void *userdata), void *userdata)
{
    self->epsilonm1 = epsilonm1;
    self->userdata  = userdata;
}

/**
 * @brief Set algorithm to calculate deterimant
 *
 * The algorithm is given by detalg. Make sure that detalg contains a valid
 * algorithm, otherwise the computation will print a warning on runtime and
 * default to DETALG_LU.
 *
 * detalg may be: DETALG_LU, DETALG_QR, DETALG_EIG
 *
 * @param [in,out] self Casimir object
 * @param [in] detalg algorithm to compute determinant
 */
void casimir_set_detalg(casimir_t *self, detalg_t detalg)
{
    self->detalg = detalg;
}

/**
 * @brief Get algorithm to calculate determinant
 *
 * The string is stored in detalg.
 *
 * This function is not thread-safe.
 *
 * @param [in,out] self Casimir object
 * @param [out] detalg buffer
 * @retval 1
 */
detalg_t casimir_get_detalg(casimir_t *self)
{
    return self->detalg;
}

/**
 * @brief Set dimension of vector space
 *
 * The round trip matrices are infinite. For a numerical evaluation the
 * dimension has to be limited to a finite value. The accuracy of the result
 * depends on the truncation of the vector space. ldim determines the dimension
 * in the angular momenta l that is used. The main contributions come from l1≈l2≈X.
 * X can be determined using \ref casimir_estimate_lminmax.
 *
 * @param [in,out] self Casimir object
 * @param [in] ldim dimension in angular momenta l
 * @retval 1 if successful
 * @retval 0 if ldim < 1
 */
int casimir_set_ldim(casimir_t *self, int ldim)
{
    if(ldim < 1)
        return 0;

    self->ldim = ldim;
    return 1;
}


/**
 * @brief Get dimension of vector space
 *
 * See \ref casimir_set_ldim.
 *
 * @param [in,out] self Casimir object
 * @retval ldim dimension of vector space
 */
int casimir_get_ldim(casimir_t *self)
{
    return self->ldim;
}


/**
 * @brief Get threshold
 *
 * See \ref casimir_set_threshold.
 *
 * @param [in,out] self Casimir object
 * @retval threshold
 */
double casimir_get_threshold(casimir_t *self)
{
    return self->threshold;
}

/**
 * @brief Set thresold
 *
 * @param [in,out] self Casimir object
 * @param [in] threshold
 * @retval 0 if threshold <= 0
 * @retval 1 if successful
 */
int casimir_set_threshold(casimir_t *self, double threshold)
{
    if(threshold < 0)
        return 0;

    self->threshold = threshold;
    return 1;
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
 * The Mie coefficients are evaluated at \f$\chi = nTR/(R+L)\f$.
 *
 * The pointers a0, sign_a0, b0 and sign_b0 must not be NULL.
 *
 * @param [in] l \f$\ell\f$
 * @param [out] a0 coefficient \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_a0 sign of \f$a_{\ell,0}^\mathrm{perf}\f$
 * @param [out] b0 coefficient \f$b_{\ell,0}^\mathrm{perf}\f$
 * @param [out] sign_b0 sign of \f$b_{\ell,0}^\mathrm{perf}\f$
 */
void casimir_lnab0(int l, double *a0, sign_t *sign_a0, double *b0, sign_t *sign_b0)
{
    *sign_a0 = MPOW(l);
    *sign_b0 = MPOW(l+1);
    *b0 = M_LOGPI-lgamma(l+0.5)-lgamma(l+1.5);
    *a0 = *b0+log1p(1.0/l);
}

/**
 * @brief Calculate Mie coefficients for perfect reflectors
 *
 * This function calculates the logarithms of the Mie coefficients
 * \f$a_\ell(i\chi)\f$ and \f$b_\ell(i\chi)\f$ for perfect reflectors and their
 * signs. The Mie coefficients are evaluated at the argument
 * \f$\chi=nT R/(R+L)\f$.
 *
 * lna, lnb, sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * @param [in,out] self Casimir object
 * @param [in] nT Matsubara frequency
 * @param [in] l angular momentum \f$\ell\f$
 * @param [out] ln_a logarithm of \f$a_\ell\f$
 * @param [out] sign sign of \f$a_\ell\f$
 * @param [out] ln_b logarithm of \f$b_\ell\f$
 * @param [out] sign sign of \f$b_\ell\f$
 */
void casimir_lnab_perf(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
{
    double lnKlp,lnKlm,lnIlm,lnIlp;
    const double chi = nT/(1+self->LbyR); /* xi*R/(R+L) = xi/(1+L/R) */

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother?
     */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp);

    /* Calculate b_l(chi), i.e. lnb and sign_b */
    *lnb = M_LOGPI-M_LOG2+lnIlp-lnKlp;

    /* We want to calculate
     * a_l(chi) = (-1)^(l+1)*pi/2 * ( l*Ip-chi*Im )/( l*Kp+chi*Km )
     *          = (-1)^(l+1)*pi/2*Ip/Kp * ( l-chi*Im/Ip )/( l+chi*Km/Kp )
     *            \--------/ \--------/   \-------------/ \-------------/
     *               sign    |b_l(chi)|      numerator      denominator
     *
     *          = b_l(chi) * numerator/denominator
     *
     * Note that chi,Km,Kp>0 and thus denominator >= 1 (and it has positive
     * sign). Also, sign_numerator is always -1 and thus the sign of al and bl
     * are always different.
     */

    /* numerator and denominator to calculate al */
    sign_t sign_numerator;
    double log_chi = log(chi);
    double numerator   = logadd_s(logi(l), +1, log_chi+lnIlm-lnIlp, -1, &sign_numerator);
    double denominator = logadd(logi(l), log_chi+lnKlm-lnKlp);

    *lna = *lnb+numerator-denominator;

    if(sign_b != NULL)
        *sign_b = MPOW(l+1);
    if(sign_a != NULL)
        *sign_a = MPOW(l+1)*sign_numerator;
}

/**
 * @brief Return logarithm of Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for arbitrary metals
 *
 * For \f$\omega_\mathrm{P} = \infty\f$ the Mie coefficient for perfect
 * reflectors are returned (see \ref casimir_lnab_perf).
 *
 * lna, lnb, sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * For generic metals we calculate the Mie coefficients al(iξ) und bl(iξ) using
 * the expressions taken from [1]. Ref. [1] is the erratum to [2]. Please note
 * that the equations (3.30) and (3.31) in [3] are wrong.
 *
 * Note: If sla =~ slb or slc =~ sld, there is a loss of significance when
 * calculating sla-slb or slc-sld.
 *
 * This function is thread safe - as long you don't change temperature, aspect
 * ratio and dielectric properties of sphere.
 *
 * The signs are given by sign_al = (-1)^l, sign_bl = (-1)^(l+1).
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
 * @param [in] nT Matsubara frequency
 * @param [in] l angular momentum \f$\ell\f$
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 * @param [out] sign_a sign of Mie coefficient \f$a_\ell\f$
 * @param [out] sign_b sign of Mie coefficient \f$b_\ell\f$
 */
void casimir_lnab(casimir_t *self, double nT, int l, double *lna, double *lnb, sign_t *sign_a, sign_t *sign_b)
{
    sign_t sign;

    /* ξ = nT */
    const double epsilonm1 = casimir_epsilonm1(self, nT);

    if(isinf(epsilonm1))
    {
        /* Mie coefficients for perfect reflectors */
        casimir_lnab_perf(self, nT, l, lna, lnb, sign_a, sign_b);
        return;
    }

    /* Mie coefficients for arbitrary metals */

    /* χ = ξ*R/(R+L) = ξ/(1+L/R) */
    const double chi    = nT/(1+self->LbyR);
    const double ln_chi = log(nT)-log1p(self->LbyR);
    const double ln_l   = logi(l);

    /**
     * Note: n is the refraction index, n_mat the Matsubara index
     * n    = sqrt(ε(ξ,ω_p,γ))
     * ln_n = ln(sqrt(ε)) = ln(ε)/2 = ln(1+(ε-1))/2 = log1p(ε-1)/2
     */
    const double ln_n = log1p(epsilonm1)/2;
    const double n    = exp(ln_n);

    double lnIlp, lnKlp, lnIlm, lnKlm, lnIlp_nchi, lnKlp_nchi, lnIlm_nchi, lnKlm_nchi;

    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp); /* I_{l+0.5}(χ), K_{l+0.5}(χ) */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm); /* K_{l-0.5}(χ), K_{l-0.5}(χ) */

    bessel_lnInuKnu(l,   n*chi, &lnIlp_nchi, &lnKlp_nchi); /* I_{l+0.5}(nχ), K_{l+0.5}(nχ) */
    bessel_lnInuKnu(l-1, n*chi, &lnIlm_nchi, &lnKlm_nchi); /* K_{l-0.5}(nχ), K_{l-0.5}(nχ) */

    double ln_gammaA, ln_gammaB, ln_gammaC, ln_gammaD;
    sign_t sign_gammaA, sign_gammaB;

    ln_gammaA = log(epsilonm1)+lnIlp_nchi + logadd_s(ln_l+lnIlp, +1, ln_chi+lnIlm, -1, &sign_gammaA);
    ln_gammaC = log(epsilonm1)+lnIlp_nchi + logadd(ln_chi+lnKlm, ln_l+lnKlp);
    ln_gammaD = ln_chi + logadd(lnIlp_nchi+lnKlm, ln_n+lnKlp+lnIlm_nchi);

    {
        double ratio_chi  = bessel_continued_fraction(l-1,chi);
        double ratio_nchi = bessel_continued_fraction(l-1,n*chi);
        double c = n/ratio_chi-1/ratio_nchi;
        ln_gammaB = ln_chi + lnIlm+lnIlm_nchi+log(fabs(c));
        sign_gammaB = copysign(1,c);
    }

    *lnb = M_LOGPI-M_LOG2 + ln_gammaB-ln_gammaD;
    *lna = M_LOGPI-M_LOG2 + logadd_s(ln_gammaA, sign_gammaA, ln_gammaB, sign_gammaB, &sign) - logadd(ln_gammaC, ln_gammaD);

    if(sign_a != NULL)
        *sign_a = sign*MPOW(l+1);
    if(sign_b != NULL)
        *sign_b = MPOW(l+1)*sign_gammaB;
}
/*@}*/


/** @brief Estimate lmin and lmax for given m and dim
 *
 * Estimate the vector space. The main contributions come from l1≈l2≈X and only
 * depends on geometry, L/R, and the quantum number m. This function calculates
 * X using the formula in the high-temperature limit and calculates lmin, lmax.
 *
 * @param [in] self Casimir object
 * @param [in] m quantum number
 * @param [out] lmin_p minimum value of l
 * @param [out] lmax_p maximum value of l
 * @retval X
 */
int casimir_estimate_lminmax(casimir_t *self, int m, size_t *lmin_p, size_t *lmax_p)
{
    const size_t dim = self->ldim;

    if(m == 0)
    {
        *lmin_p = 1;
        *lmax_p = 1+dim;
        return 0;
    }

    int l;
    const double x = 1/(1+self->LbyR);

    /* find maximum, i.e., main contributions */
    double last = 2*m*log(x/2);
    for(l = m+1; 1; l++)
    {
        double f = lfac(2*l)-lfac(l+m)-lfac(l-m)+2*l*log(x/2);
        if(f < last)
        {
            l--;
            break;
        }

        last = f;
    }

    int lmin = l-dim/2;
    if(lmin < m)
        lmin = m;

    *lmin_p = lmin;
    *lmax_p = lmin + dim;

    return l;
}

/**
 * @brief Calculate matrix elements for xi=0
 *
 * This function calculates the matrix elements for <l1,m,E|M|l2,m,E> and
 * <l1,m,M|M|l2,m,M> in the high-temperature limit, i.e. xi=0.
 *
 * If EE or MM is NULL, it is not written. Otherwise the result is written to
 * memory pointed by EE and MM.
 *
 * @param [in]  self Casimir object
 * @param [in]  l1
 * @param [in]  l2
 * @param [in]  m
 * @param [out] EE
 * @param [out] MM
 */
void casimir_M0_elem(casimir_t *self, int l1, int l2, int m, double *EE, double *MM)
{
    /* y = log(R/(R+L)/2) */
    double y = self->y;

    /* See thesis of Antoine, section 6.7:
     * x = R/(R+L)
     * M_EE_{l1,l2} = (x/2)^(l1+l2) * (l1+l2)! / sqrt( (l1+m)!*(l1-m)! * (l2+m)!*(l2-m)! )
     * M_MM_{l1,l2} = M_EE_{l1,l2} * sqrt( l1/(l1+1) * l2/(l2+1) )
     */
    double _EE = exp( (l1+l2+1)*y + lfac(l1+l2) - 0.5*(lfac(l1+m)+lfac(l1-m) + lfac(l2+m)+lfac(l2-m)) );
    if(EE != NULL)
        *EE = _EE;
    if(MM != NULL)
        *MM = _EE*sqrt((l1*l2)/((l1+1.)*(l2+1.)));
}

double casimir_kernel_M0_EE(int i, int j, void *args_)
{
    casimir_M_t *args = (casimir_M_t *)args_;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;
    casimir_t *casimir = args->casimir;

    double EE;
    casimir_M0_elem(casimir, l1, l2, m, &EE, NULL);
    return EE;
}

double casimir_kernel_M0_MM(int i, int j, void *args_)
{
    casimir_M_t *args = (casimir_M_t *)args_;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;
    casimir_t *casimir = args->casimir;

    double MM;
    casimir_M0_elem(casimir, l1, l2, m, NULL, &MM);
    return MM;
}

/**
 * @brief Calculate round-trip matrices M for xi=nT=0
 *
 * For xi=0 the round-trip matrix M is block diagonal with block matrices EE,
 * MM. This function calculates the block matrices EE and MM.
 *
 * If EE is not NULL, create and calculate block matrix EE.
 * If MM is not NULL, create and calculate block matrix MM.

 * You have to free the matrices yourself using matrix_free.
 *
 * @param [in] self Casimir object
 * @param [out] EE pointer for block matrix EE
 * @param [out] MM pointer for block matrix MM
 */
void casimir_M0(casimir_t *self, int m, matrix_t **EE, matrix_t **MM)
{
    /* main contributions comes from l1≈l2≈m */
    size_t lmin,lmax,ldim = self->ldim;
    casimir_estimate_lminmax(self, m, &lmin, &lmax);

    /* nothing to do... */
    if(EE == NULL && MM == NULL)
        return;

    if(EE != NULL)
    {
        *EE = matrix_alloc(ldim);
        matrix_setall(*EE,0);
    }
    if(MM != NULL)
    {
        *MM = matrix_alloc(ldim);
        matrix_setall(*MM,0);
    }

    double trace_diag = 0; /* sum of modulus of matrix elements of diagonal */

    /* calculate matrix elements of M */

    /* n-th minor diagonal */
    for(size_t n = 0; n < ldim; n++)
    {
        /* sum of modulus of matrix elements of n-th minor diagonal */
        double trace = 0;

        for(size_t d = 0; d < ldim-n; d++)
        {
            double EE_ij, MM_ij;
            const int l1 = d+lmin;
            const int l2 = d+n+lmin;

            /* i: row of matrix, j: column of matrix */
            const size_t i = d, j = d+n;

            casimir_M0_elem(self, l1, l2, m, &EE_ij, &MM_ij);

            /* calculate trace of n-th minor diagonal */
            trace += fabs(EE_ij);

            /* The matrix M is symmetric. */
            if(EE != NULL)
            {
                matrix_set(*EE, i,j, EE_ij);
                matrix_set(*EE, j,i, EE_ij);
            }
            if(MM != NULL)
            {
                matrix_set(*MM, i,j, MM_ij);
                matrix_set(*MM, j,i, MM_ij);
            }
        }

        if(n == 0)
            trace_diag = ldim-trace;
        else if(trace/trace_diag <= self->threshold)
            break;
    }
}

/**
 * @brief Calculate \f$\log\det \mathcal{D}^{(m)}(\xi=0)\f$ for EE and MM
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n=0\f$.
 *
 * This function is thread-safe as long you don't change ldim and aspect ratio.
 *
 * @param [in,out] self Casimir object
 * @param [in] m
 * @param [out] logdet_EE
 * @param [out] logdet_MM
 */
void casimir_logdetD0_dense(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    matrix_t *EE = NULL, *MM = NULL;

    if(logdet_EE != NULL && logdet_MM != NULL)
        casimir_M0(self, m, &EE, &MM);
    else if(logdet_EE == NULL && logdet_MM != NULL)
        casimir_M0(self, m, NULL, &MM);
    else if(logdet_EE != NULL && logdet_MM == NULL)
        casimir_M0(self, m, &EE, NULL);

    /* Calculate logdet=log(det(Id-M)) and free space. */
    if(EE != NULL)
    {
        *logdet_EE = matrix_logdet(EE, -1);
        matrix_free(EE);
    }
    if(MM != NULL)
    {
        *logdet_MM = matrix_logdet(MM, -1);
        matrix_free(MM);
    }
}

void casimir_logdetD0(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    if(self->detalg == DETALG_HODLR)
        casimir_logdetD0_hodlr(self, m, logdet_EE, logdet_MM);
    else
        casimir_logdetD0_dense(self, m, logdet_EE, logdet_MM);
}

void casimir_logdetD0_hodlr(casimir_t *self, int m, double *logdet_EE, double *logdet_MM)
{
    size_t lmin, lmax;
    unsigned int nLeaf = 100; /* XXX */
    int is_symmetric = 1;
    double tolerance = 1e-15;

    casimir_estimate_lminmax(self, m, &lmin, &lmax);

    casimir_M_t args = {
        .casimir = self,
        .m = m,
        .nT = 0,
        .integration = NULL,
        .al = NULL,
        .bl = NULL,
        .lmin = lmin
    };

    if(logdet_EE != NULL)
        *logdet_EE = hodlr_logdet(self->ldim, &casimir_kernel_M0_EE, &args, nLeaf, tolerance, is_symmetric);

    if(logdet_MM != NULL)
        *logdet_MM = hodlr_logdet(self->ldim, &casimir_kernel_M0_MM, &args, nLeaf, tolerance, is_symmetric);
}

double casimir_kernel_M(int i, int j, void *args_)
{
    char p1 = 'E', p2 = 'E';
    casimir_M_t *args = (casimir_M_t *)args_;
    const int ldim = args->casimir->ldim;
    const int lmin = args->lmin;
    int l1 = i+lmin, l2 = j+lmin;

    if(i >= ldim)
    {
        p1 = 'M';
        l1 -= ldim;
    }
    if(j >= ldim)
    {
        p2 = 'M';
        l2 -= ldim;
    }

    return casimir_M_elem(args, l1, l2, p1, p2);
}

casimir_M_t *casimir_M_init(casimir_t *casimir, int m, double nT)
{
    size_t lmin, lmax;
    const int ldim = casimir->ldim;
    casimir_M_t *self = xmalloc(sizeof(casimir_M_t));

    casimir_estimate_lminmax(casimir, m, &lmin, &lmax);

    self->casimir = casimir;
    self->m = m;
    self->lmin = lmin; // XXX
    self->integration = casimir_integrate_init(casimir, nT, m, casimir->tolerance);
    self->nT = nT;
    self->al = xmalloc(ldim*sizeof(double));
    self->bl = xmalloc(ldim*sizeof(double));

    for(int j = 0; j < ldim; j++)
        self->al[j] = self->bl[j] = NAN;

    return self;
}

double casimir_M_elem(casimir_M_t *self, int l1, int l2, char p1, char p2)
{
    const double nT = self->nT;
    const int lmin = self->lmin;
    casimir_t *casimir = self->casimir;
    integration_t *integration = self->integration;

    if(isnan(self->al[l1-lmin]))
        casimir_lnab(casimir, nT, l1, &self->al[l1-lmin], &self->bl[l1-lmin], NULL, NULL);

    if(isnan(self->al[l2-lmin]))
        casimir_lnab(casimir, nT, l2, &self->al[l2-lmin], &self->bl[l2-lmin], NULL, NULL);

    const double al1 = self->al[l1-lmin], bl1 = self->bl[l1-lmin];
    const double al2 = self->al[l2-lmin], bl2 = self->bl[l2-lmin];

    if(p1 == p2) /* EE or MM */
    {
        if(p1 == 'E') /* EE */
        {
            sign_t signA_TE, signB_TM;
            double log_A_TE = casimir_integrate_A(integration, l1, l2, TE, &signA_TE);
            double log_B_TM = casimir_integrate_B(integration, l1, l2, TM, &signB_TM);

            /* √(a_l1*a_l2)*(A_TE + B_TM) */
            const double mie = (al1+al2)/2;
            const double elem = exp(log_A_TE+mie)*signA_TE+exp(log_B_TM+mie)*signB_TM;

            return MPOW(l1)*elem;
        }
        else /* MM */
        {
            sign_t signA_TM, signB_TE;
            double log_A_TM = casimir_integrate_A(integration, l1, l2, TM, &signA_TM);
            double log_B_TE = casimir_integrate_B(integration, l1, l2, TE, &signB_TE);

            /* √(b_l1*b_l2)*(A_TM + B_TE) */
            const double mie = (bl1+bl2)/2;
            const double elem = exp(log_A_TM+mie)*signA_TM+exp(log_B_TE+mie)*signB_TE;

            return MPOW(l1+1)*elem;
        }
    }
    else /* EM or ME */
    {
        const int m = self->m;
        if(m == 0)
            return 0;

        /* M_EM */
        if(p1 == 'E') /* EM */
        {
            sign_t signC_TE, signD_TM;

            double log_C_TE = casimir_integrate_C(integration, l1, l2, TE, &signC_TE);
            double log_D_TM = casimir_integrate_D(integration, l1, l2, TM, &signD_TM);

            /* C_TE + D_TM */
            double mie1 = (al1+bl2)/2;
            double elem1 = exp(log_C_TE+mie1)*signC_TE+exp(log_D_TM+mie1)*signD_TM;

            return MPOW(l1)*elem1;
        }
        else /* ME */
        {
            sign_t signC_TM, signD_TE;

            double log_C_TM = casimir_integrate_C(integration, l1, l2, TM, &signC_TM);
            double log_D_TE = casimir_integrate_D(integration, l1, l2, TE, &signD_TE);

            /* C_TM + D_TE */
            const double mie2 = (bl1+al2)/2;
            const double elem2 = exp(log_C_TM+mie2)*signC_TM+exp(log_D_TE+mie2)*signD_TE;

            return -MPOW(l1)*elem2;
        }
    }
}

void casimir_M_free(casimir_M_t *self)
{
    xfree(self->al);
    xfree(self->bl);
    casimir_integrate_free(self->integration);
    xfree(self);
}

/**
 * @brief Calculate round-trip matrix M
 *
 * Create the round-trip matrix M for Matsubara frequency xi=nT and angular
 * momentum number m.
 *
 * You have to free the matrix yourself using matrix_free.
 *
 * @param [in] self Casimir object
 * @param [in] nT Matsubara frequency
 * @param [in] m
 * @retval M round-trip matrix
 */
matrix_t *casimir_M(casimir_t *self, double nT, int m)
{
    /* main contributions comes from l1≈l2 */
    size_t lmin, lmax, ldim = self->ldim;
    casimir_estimate_lminmax(self, m, &lmin, &lmax);

    casimir_M_t *obj = casimir_M_init(self, m, nT);

    /* allocate space for matrix M */
    matrix_t *M = matrix_alloc(2*ldim);

    /* set matrix elements to 0 */
    matrix_setall(M, 0);

    /* M_EE, -M_EM
       M_ME,  M_MM */

    double trace_diag = 0;

    /* n-th minor diagonal */
    for(size_t md = 0; md < ldim; md++)
    {
        double trace = 0;

        for(size_t k = 0; k < ldim-md; k++)
        {
            const int l1 = k+lmin;
            const int l2 = k+md+lmin;

            /* i: row of matrix, j: column of matrix */
            const size_t i = l1-lmin, j = l2-lmin;

            /* EE */
            double EE = casimir_M_elem(obj, l1, l2, 'E', 'E');
            matrix_set(M, i,j, EE);
            matrix_set(M, j,i, EE);

            /* MM */
            double MM = casimir_M_elem(obj, l1, l2, 'M', 'M');
            matrix_set(M, i+ldim,j+ldim, MM);
            matrix_set(M, j+ldim,i+ldim, MM);

            trace += fabs(EE)+fabs(MM);

            if(m)
            {
                /* non-diagonal blocks EM and ME */
                double EM = casimir_M_elem(obj, l1, l2, 'E', 'M');
                double ME = casimir_M_elem(obj, l1, l2, 'M', 'E');

                /* EM */
                matrix_set(M, i,ldim+j, EM);
                matrix_set(M, j,ldim+i, ME);

                /* ME */
                matrix_set(M, ldim+i,j, ME);
                matrix_set(M, ldim+j,i, EM);
            }

            #if 0
            /* ME */
            if(m != 0)
            {
                sign_t signC_TE, signC_TM, signD_TE, signD_TM;

                const double log_C_TE = casimir_integrate_C(integration, l1, l2, TE, &signC_TE);
                const double log_C_TM = casimir_integrate_C(integration, l1, l2, TM, &signC_TM);

                const double log_D_TE = casimir_integrate_D(integration, l1, l2, TE, &signD_TE);
                const double log_D_TM = casimir_integrate_D(integration, l1, l2, TM, &signD_TM);

                const double mie1 = (ln_al1+ln_bl2)/2;
                const double mie2 = (ln_bl1+ln_al2)/2;

                /* C_TE + D_TM */
                const double elem1 = exp(log_C_TE+mie1)*signC_TE+exp(log_D_TM+mie1)*signD_TM;

                /* C_TM + D_TE */
                const double elem2 = exp(log_C_TM+mie2)*signC_TM+exp(log_D_TE+mie2)*signD_TE;

                /* M_EM */
                matrix_set(M, i,ldim+j,  MPOW(l1)*elem1);
                matrix_set(M, j,ldim+i, -MPOW(l1)*elem2);

                /* M_ME */
                matrix_set(M, ldim+i,j, -MPOW(l1)*elem2);
                matrix_set(M, ldim+j,i,  MPOW(l1)*elem1);
            }
            #endif
        }

        if(md == 0)
            trace_diag = 2*ldim - trace;
        else if(trace/trace_diag <= self->threshold)
            break;
    }

    casimir_M_free(obj);

    return M;
}

double casimir_logdetD(casimir_t *self, double nT, int m)
{
    if(self->detalg == DETALG_HODLR)
        return casimir_logdetD_hodlr(self, nT, m);
    else
        return casimir_logdetD_dense(self, nT, m);
}

double casimir_logdetD_hodlr(casimir_t *self, double nT, int m)
{
    unsigned int nLeaf = 100; /* XXX */
    int is_symmetric = 1;
    double tolerance = 1e-15;

    casimir_M_t *obj = casimir_M_init(self, m, nT);
    double logdet = hodlr_logdet(2*self->ldim, &casimir_kernel_M, obj, nLeaf, tolerance, is_symmetric);

    casimir_M_free(obj);

    return logdet;
}

/**
 * @brief Calculate \f$\log\det \mathcal{D}^{(m)}(\xi=nT)\f$
 *
 * This function calculates the logarithm of the determinant of the scattering
 * operator D for the Matsubara term \f$n\f$.
 *
 * This function is thread-safe - as long you don't change ldim, temperature,
 * aspect ratio, dielectric properties of sphere and plane, and integration.
 *
 * @param [in,out] self Casimir object
 * @param [in] nT Matsubara frequency
 * @param [in] m
 * @retval logdetD \f$\log \det \mathcal{D}^{(m)}(\xi=nT)\f$
 */
double casimir_logdetD_dense(casimir_t *self, double nT, int m)
{
    double t0;
    double logdet = 0;

    if(nT == 0)
    {
        double logdet_EE = 0, logdet_MM = 0;

        /* XXX what happens for xi=0 XXX */
        //if(self->epsilonm1 == NULL)
            /* perfect reflector */
            casimir_logdetD0(self, m, &logdet_EE, &logdet_MM);
        //else
        //    /* generic mirrors */
        //    casimir_logdetD0(self, m, &logdet_EE, NULL);

        return logdet_EE + logdet_MM;
    }

    t0 = now();
    matrix_t *M = casimir_M(self, nT, m);
    casimir_debug(self, "# timing: matrix elements: %gs\n", now()-t0);

    #if 0
    /* dump matrix */
    t0 = now();
    matrix_save_to_file(M, "M.npy");
    casimir_debug(self, "# dumped round-trip matrix M %g\n", now()-t0);
    #endif

    if(m == 0)
    {
        const size_t dim  = M->dim;
        const size_t ldim = dim/2;

        /* create a view of the block matrices EE and MM */
        matrix_t *EE = matrix_view(M->M, ldim, dim);
        matrix_t *MM = matrix_view(&M->M[ldim*(dim+1)], ldim, dim);

        t0 = now();
        logdet  = matrix_logdet(EE, -1);
        logdet += matrix_logdet(MM, -1);
        casimir_debug(self, "# timing: log(det(Id-EE)), log(det(Id-MM)): %gs\n", now()-t0);

        matrix_free(EE);
        matrix_free(MM);
    }
    else
    {
        t0 = now();
        logdet = matrix_logdet(M, -1);
        casimir_debug(self, "# timing: log(det(Id-M)): %gs\n", now()-t0);
    }

    matrix_free(M);

    return logdet;
}

/*@}*/
