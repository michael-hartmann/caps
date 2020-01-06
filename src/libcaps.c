/**
 * @file   libcaps.c
 * @author Michael Hartmann <caps@speicherleck.de>
 * @date   December, 2017
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "quadpack.h"
#include "bessel.h"
#include "integration.h"
#include "libcaps.h"
#include "matrix.h"
#include "logfac.h"
#include "misc.h"
#include "utils.h"

/**
* @name various functions
*/
/*@{*/

/** @brief Calculate logarithm \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$ for
 * \f$\ell_1,\ell_2,m\f$.
 * \f[
 *      \Lambda_{\ell_1,\ell_2}^{(m)} = \frac{2 N_{\ell_1,m} N_{\ell_2,m}}{\sqrt{\ell_1 (\ell_1+1) \ell_2 (\ell_2+1)}}
 * \f]
 *
 * Symmetries: \f$\Lambda_{\ell_1,\ell_2}^{(m)} = \Lambda_{\ell_2,\ell_1}^{(m)}\f$
 *
 * @param [in]  l1 l1>0
 * @param [in]  l2 l2>0
 * @param [in]  m  m <= l1 and m <= l2
 * @retval lnLambda \f$\log{\Lambda_{\ell_1,\ell_2}^{(m)}}\f$
 */
double caps_lnLambda(int l1, int l2, int m)
{
    return (logi(2*l1+1)+logi(2*l2+1)-logi(l1)-logi(l1+1)-logi(l2)-logi(l2+1)+lfac(l1-m)+lfac(l2-m)-lfac(l1+m)-lfac(l2+m))/2.0;
}

/** @brief Estimate \f$\ell_\mathrm{min}\f$ and \f$\ell_\mathrm{max}\f$
 *
 * Estimate the vector space: The main contributions comes from the vicinity
 * \f$\ell_1=\ell_2=X\f$ and only depend on geometry, \f$L/R\f$, and the quantum number \f$m\f$. This
 * function calculates \f$X\f$ using the formula in the high-temperature limit and
 * calculates \f$\ell_\mathrm{min}\f$, \f$\ell_\mathrm{max}\f$.
 *
 * @param [in] self CaPS object
 * @param [in] m quantum number
 * @param [out] lmin_p minimum value of \f$\ell\f$
 * @param [out] lmax_p maximum value of \f$\ell\f$
 * @retval l approximately the value of \f$\ell\f$ where \f$\mathcal{M}^m_{\ell\ell}\f$ is maximal
 */
int caps_estimate_lminmax(caps_t *self, int m, size_t *lmin_p, size_t *lmax_p)
{
    const int ldim = self->ldim;

    if(m == 0)
    {
        *lmin_p = 1;
        *lmax_p = 1+ldim;
        return 0;
    }

    /* approximately the maximum on the diagonal of the round-trip matrix M is
     * at l */
    const int l = round(m/sqrt(2*self->LbyR));

    /* lmin >= m */
    int lmin = l-ldim/2;
    if(lmin < m)
        lmin = m;

    /* minimum and maximum l */
    *lmin_p = lmin;
    *lmax_p = lmin + ldim;

    return l;
}


/*@}*/

/**
* @name Dielectric functions
*/
/*@{*/

/** @brief Evaluate dielectric function of the plate
 *
 * @param [in] self CaPS object
 * @param [in] xi_ \f$\xi\mathcal{L}/c\f$
 * @retval epsm1 \f$\epsilon(\mathrm{i}\xi)\f$
 */
double caps_epsilonm1_plate(caps_t *self, double xi_)
{
    double xi = xi_*CAPS_C/self->calL;
    return self->epsilonm1_plate(xi, self->userdata_plate);
}

/** @brief Evaluate dielectric function of the sphere
 *
 * @param [in] self CaPS object
 * @param [in] xi_ \f$\xi\mathcal{L}/c\f$
 * @retval epsm1 \f$\epsilon(\mathrm{i}\xi)\f$
 */
double caps_epsilonm1_sphere(caps_t *self, double xi_)
{
    double xi = xi_*CAPS_C/self->calL;
    return self->epsilonm1_sphere(xi, self->userdata_sphere);
}

/**
 * @brief Dielectric function for perfect reflectors
 *
 * @param [in] xi_ ignored
 * @param [in] userdata ignored
 * @retval inf \f$\epsilon(\xi) = \infty\f$
 */
double caps_epsilonm1_perf(__attribute__((unused)) double xi_, __attribute__((unused)) void *userdata)
{
    return INFINITY;
}

/**
 * @brief Dielectric function for Drude reflectors
 *
 * Dielectric function for Drude
 * \f[
 *      \epsilon(\xi)-1 = \frac{\omega_\mathrm{P}^2}{\xi(\xi+\gamma)}
 * \f]
 *
 * The parameters \f$\omega_\mathrm{P}\f$ and \f$\gamma\f$ must be provided by userdata:
 *  - userdata[0] = \f$\omega_\mathrm{P}\f$ in rad/s
 *  - userdata[1] = \f$\gamma\f$ in rad/s
 *
 * @param [in] xi frequency in rad/s
 * @param [in] userdata userdata
 * @retval epsilon epsilon(xi)
 */
double caps_epsilonm1_drude(double xi, void *userdata)
{
    double *ptr = (double *)userdata;
    double omegap = ptr[0], gamma_ = ptr[1];

    return POW_2(omegap)/(xi*(xi+gamma_));
}

/*@}*/

/**
* @name initialization and setting parameters
*/
/*@{*/

/**
 * @brief Create a new CaPS object
 *
 * This function will initialize a CaPS object. By default the dielectric
 * function corresponds to perfect reflectors, i.e. \f$\epsilon(\xi)=\infty\f$.
 *
 * By default, the value of \f$\ell_\mathrm{dim}\f$ is chosen by:
 * \f[
 *  \ell_\mathrm{dim} = \mathrm{ceil}\left(\mathrm{max}\left(\mathrm{CAPS\_MINIMUM\_LDIM}, \mathrm{CAPS\_FACTOR\_LDIM}\cdot \frac{R}{L}\right)\right)
 * \f]
 *
 * Restrictions: \f$L/R > 0\f$
 *
 * @param [in]  R radius of sphere in m
 * @param [in]  L smallest separation between sphere and plate in m
 * @retval object CaPS object if successful
 * @retval NULL   if an error occured
 */
caps_t *caps_init(double R, double L)
{
    const double LbyR = L/R;
    if(L <= 0 || R <= 0)
        return NULL;

    caps_t *self = xmalloc(sizeof(caps_t));

    /* geometry */
    self->LbyR = LbyR;
    self->L    = L;
    self->R    = R;
    self->calL = L+R;
    self->y = -CAPS_LOG2-log1p(LbyR); /* log( (R/(L+R))/2 ) */

    /* dimension of vector space */
    self->ldim = ceil(MAX(CAPS_MINIMUM_LDIM, CAPS_FACTOR_LDIM/LbyR));

    /* perfect reflectors */
    self->epsilonm1_plate  = caps_epsilonm1_perf;
    self->epsilonm1_sphere = caps_epsilonm1_perf;
    self->userdata_plate   = NULL;
    self->userdata_sphere  = NULL;

    /* relative error for integration */
    self->epsrel = CAPS_EPSREL;

    /* use LU decomposition by default */
    self->detalg = DETALG_HODLR;

    return self;
}


/**
 * @brief Free memory for CaPS object
 *
 * Free allocated memory for the CaPS object self.
 *
 * @param [in,out] self CaPS object
 */
void caps_free(caps_t *self)
{
    xfree(self);
}

/**
 * @brief Print information on build to stream
 *
 * The information contains compiler, build time, git head and git branch if
 * available. If prefix is not NULL, the string prefix will added in front of
 * each line.
 *
 * @param stream output stream
 * @param prefix prefix of each line or NULL
 */
void caps_build(FILE *stream, const char *prefix)
{
    if(prefix == NULL)
        prefix = "";

    #ifdef LIBCAPS_VERSION
    fprintf(stream, "%sversion: %s\n", prefix, LIBCAPS_VERSION);
    #endif

    fprintf(stream, "%scompiler: %s\n", prefix, COMPILER);
    fprintf(stream, "%scompile time: %s %s\n", prefix, __DATE__, __TIME__);

    #ifdef MACHINE
    fprintf(stream, "%scompiled on: %s\n", prefix, MACHINE);
    #endif

    #ifdef GIT_HEAD
    fprintf(stream, "%sgit HEAD: %s\n", prefix, GIT_HEAD);
    #endif

    #ifdef GIT_BRANCH
    fprintf(stream, "%sgit branch: %s\n", prefix, GIT_BRANCH);
    #endif
}

/** @brief Print object information to stream
 *
 * Print information about the object self to stream.
 *
 * @param self CaPS object
 * @param stream where to print the string
 * @param prefix if prefix != NULL: start every line with the string contained
 * in prefix
 */
void caps_info(caps_t *self, FILE *stream, const char *prefix)
{
    const char *detalg_str;
    if(prefix == NULL)
        prefix = "";

    caps_detalg_to_string(self->detalg, &detalg_str);

    fprintf(stream, "%sL/R    = %.16g\n", prefix, self->LbyR);
    fprintf(stream, "%sL      = %.16g\n", prefix, self->L);
    fprintf(stream, "%sR      = %.16g\n", prefix, self->R);
    fprintf(stream, "%sldim   = %d\n",    prefix, self->ldim);
    fprintf(stream, "%sepsrel = %.1e\n",  prefix, self->epsrel);
    fprintf(stream, "%sdetalg = %s\n",    prefix, detalg_str);
}

/**
 * @brief Set relative error for numerical integration
 *
 * Set relative error for numerical integration.
 *
 * @param [in] self CaPS object
 * @param [in] epsrel relative error
 * @retval 0 if an error occured
 * @retval 1 on success
 */
int caps_set_epsrel(caps_t *self, double epsrel)
{
    if(epsrel <= 0)
        return 0;

    self->epsrel = epsrel;
    return 1;
}

/**
 * @brief Get relative error for numerical integration
 *
 * See \ref caps_set_epsrel.
 *
 * @retval epsrel relative error
 */
double caps_get_epsrel(caps_t *self)
{
    return self->epsrel;
}

/**
 * @brief Set dielectric function for plate and sphere
 *
 * See also \ref caps_set_epsilonm1_plate and \ref
 * caps_set_epsilonm1_sphere.
 *
 * @param [in,out] self CaPS object
 * @param [in] epsilonm1  callback to the function that calculates \f$\epsilon(\mathrm{i}\xi)-1\f$
 * @param [in] userdata   arbitrary pointer to data that is passwd to epsilonm1 whenever the function is called
 */
void caps_set_epsilonm1(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata)
{
    caps_set_epsilonm1_plate (self, epsilonm1, userdata);
    caps_set_epsilonm1_sphere(self, epsilonm1, userdata);
}

/**
 * @brief Set dielectric function of plate
 *
 * The Fresnel coefficient \f$r_p\f$ depend on the dielectric
 * function \f$\epsilon(\mathrm{i}\xi)\f$. By default, perfect reflectors with
 * a dielectric function \f$\epsilon(\mathrm{i}\xi)=\infty\f$ are used.
 *
 * However, you can also specify an arbitrary function for
 * \f$\epsilon(\mathrm{i}\xi)\f$. userdata is an arbitrary pointer that will be
 * given to the callback function.
 *
 * @param [in,out] self CaPS object
 * @param [in] epsilonm1  callback to the function that calculates \f$\epsilon(\mathrm{i}\xi)-1\f$
 * @param [in] userdata   arbitrary pointer to data that is passwd to epsilonm1 whenever the function is called
 */
void caps_set_epsilonm1_plate(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata)
{
    self->epsilonm1_plate = epsilonm1;
    self->userdata_plate  = userdata;
}

/**
 * @brief Set dielectric function of sphere
 *
 * The Mie coefficient \f$a_\ell,b_\ell\f$ depend on the dielectric function
 * \f$\epsilon(\mathrm{i}\xi)\f$. By default, perfect reflectors with a
 * dielectric function \f$\epsilon(\mathrm{i}\xi)=\infty\f$ are used.
 *
 * However, you can also specify an arbitrary function for
 * \f$\epsilon(\mathrm{i}\xi)\f$. userdata is an arbitrary pointer that will be
 * given to the callback function.
 *
 * @param [in,out] self CaPS object
 * @param [in] epsilonm1  callback to the function that calculates \f$\epsilon(\mathrm{i}\xi)-1\f$
 * @param [in] userdata   arbitrary pointer to data that is passwd to epsilonm1 whenever the function is called
 */
void caps_set_epsilonm1_sphere(caps_t *self, double (*epsilonm1)(double xi_, void *userdata), void *userdata)
{
    self->epsilonm1_sphere = epsilonm1;
    self->userdata_sphere  = userdata;
}

/**
 * @brief Get detalg from string
 *
 * Given the string str, find the corresponding detalg_t value and write it to
 * detalg.
 *
 * If detalg is NULL, no value will be written to detalg.
 *
 * @param [in] str    string
 * @param [in] detalg determinant algorithm
 * @retval 1          if determinant algorithm was found
 * @retval 0          if determinant algorithm was not found
 */
int caps_detalg_from_string(const char *str, detalg_t *detalg)
{
    if(str == NULL)
        return 0;

    if(strcaseequal(str, "HODLR"))
    {
        if(detalg != NULL)
            *detalg = DETALG_HODLR;
        return 1;
    }
    else if(strcaseequal(str, "Cholesky"))
    {
        if(detalg != NULL)
            *detalg = DETALG_CHOLESKY;
        return 1;
    }
    else if(strcaseequal(str, "LU"))
    {
        if(detalg != NULL)
            *detalg = DETALG_LU;
        return 1;
    }
    else if(strcaseequal(str, "QR"))
    {
        if(detalg != NULL)
            *detalg = DETALG_QR;
        return 1;
    }

    return 0;
}

/**
 * @brief Get string description of detalg
 *
 * Given a detalg, return a string description of the algorithm.
 *
 * @param [in] detalg determinant algorithm
 * @param [in] str    pointer to a constant string
 * @retval 1          if determinant algorithm was found
 * @retval 0          if determinant algorithm was not found
 */
int caps_detalg_to_string(detalg_t detalg, const char **str)
{
    static const char *detalg_str[] = {
        "HODLR", "Cholesky", "LU", "QR", "unknown"
    };

    switch(detalg)
    {
        case DETALG_HODLR:    *str = detalg_str[0]; return 1;
        case DETALG_CHOLESKY: *str = detalg_str[1]; return 1;
        case DETALG_LU:       *str = detalg_str[2]; return 1;
        case DETALG_QR:       *str = detalg_str[3]; return 1;
        default:              *str = detalg_str[4]; return 0;
    }
}

/**
 * @brief Set algorithm to calculate determinant
 *
 * The algorithm is given by detalg. Usually you don't want to change the
 * algorithm to compute the determinant.
 *
 * detalg may be: DETALG_HODLR or DETALG_LU, DETALG_QR, DETALG_CHOLESKY.
 *
 * If successul, the function returns 1. If the algorithm is not supported
 * because of missing LAPACK support, 0 is returned.
 *
 * @param [in,out] self CaPS object
 * @param [in] detalg algorithm to compute determinant
 */
void caps_set_detalg(caps_t *self, detalg_t detalg)
{
    self->detalg = detalg;
}

/**
 * @brief Get algorithm to calculate determinant
 *
 * @param [in]  self CaPS object
 * @retval detalg
 */
detalg_t caps_get_detalg(caps_t *self)
{
    return self->detalg;
}

/**
 * @brief Set dimension of vector space
 *
 * The round trip matrices are infinite. For a numerical evaluation the
 * dimension has to be truncated to a finite value. The accuracy of the result
 * depends on the truncation of the vector space. ldim determines the dimension
 * in the angular momentum \f$\ell\f$ that is used.
 *
 * @param [in,out] self CaPS object
 * @param [in] ldim dimension in angular momentum \f$\ell\f$
 * @retval 1 if successful
 * @retval 0 if ldim < 1
 */
int caps_set_ldim(caps_t *self, int ldim)
{
    if(ldim < 1)
        return 0;

    self->ldim = ldim;
    return 1;
}


/**
 * @brief Get dimension of vector space
 *
 * See \ref caps_set_ldim.
 *
 * @param [in,out] self CaPS object
 * @retval ldim dimension of vector space
 */
int caps_get_ldim(caps_t *self)
{
    return self->ldim;
}

/*@}*/


/**
 * @name Mie and Fresnell coefficients
 */
/*@{*/

/**
 * @brief Calculate Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for perfect reflectors
 *
 * This function calculates the logarithms of the Mie coefficients
 * \f$a_\ell(i\chi)\f$ and \f$b_\ell(i\chi)\f$ for perfect reflectors. The Mie
 * coefficients are evaluated at the argument \f$\chi=\xi R/c\f$.
 *
 * The signs are given by \f$\mathrm{sgn}(a_\ell) = (-1)^\ell\f$,
 * \f$\mathrm{sgn}(b_\ell) = (-1)^{\ell+1}\f$.
 *
 * lna and lnb must be valid pointers and must not be NULL.
 *
 * @param [in,out] self CaPS object
 * @param [in] xi_ \f$\xi\mathcal{L}/c > 0\f$
 * @param [in] l angular momentum \f$\ell > 0\f$
 * @param [out] lna logarithm of \f$|a_\ell|\f$
 * @param [out] lnb logarithm of \f$|b_\ell|\f$
 */
void caps_mie_perf(caps_t *self, double xi_, int l, double *lna, double *lnb)
{
    double logKlp,logKlm,logIlm,logIlp;
    /* χ = ξR/c = ξ(R+L)/c * R/(R+L) = xi_ 1/(1+L/R) */
    const double chi = xi_/(1+self->LbyR);
    const double log_chi = log(xi_)-log1p(self->LbyR);

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother?
     */
    bessel_logInKn_half(l-1, chi, &logIlm, &logKlm);
    bessel_logInKn_half(l,   chi, &logIlp, &logKlp);

    /* Calculate b_l(chi), i.e. lnb */
    *lnb = CAPS_LOGPI-CAPS_LOG2+logIlp-logKlp;

    /* We want to calculate
     *
     * b_l(χ) = π/2 * Ip/Im
     * a_l(χ) = π/2 * ( χ*Ilm - l*Ilp )/( l*Kp + χ*Km )
     *        = π/2 * Ilp * ( χ*Ilm/Ilp - l )/( l*Kp + χ*Km )
     *                          \-----/
     *                           ratio
     *
     * where Ip = I_{l+1/2}(χ), Im = I_{l-1/2}(χ), and similar for Kp and Km.
     *
     * Also note that all terms in brackets are positive.
     */

    /* numerator and denominator to calculate al */
    double ratio = bessel_ratioI(l-0.5, chi);
    double numerator = CAPS_LOGPI-CAPS_LOG2 + logIlp + log(chi*ratio - l);
    double denominator = logadd(logi(l)+logKlp, log_chi+logKlm);

    *lna = numerator-denominator;
}

/**
 * @brief Return logarithm of Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for arbitrary metals
 *
 * For \f$\omega_\mathrm{P} = \infty\f$ the Mie coefficient for perfect
 * reflectors are returned (see \ref caps_mie_perf).
 *
 * lna and lnb must be valid pointers.
 *
 * For generic metals, we calculate the Mie coefficients \f$a_\ell\f$ und
 * \f$b_\ell\f$ using the expressions taken from [1]. Ref. [1] is the erratum
 * to [2]. Please note that the equations (3.30) and (3.31) in [3] are wrong.
 * The formulas are corrected in [4].
 *
 * Note: If sla \f$\approx\f$ slb or slc \f$\approx\f$ sld, there is a loss of
 * significance when calculating sla-slb or slc-sld.
 *
 * The signs are given by \f$\mathrm{sgn}(a_\ell) = (-1)^\ell\f$, \f$\mathrm{sgn}(b_\ell) = (-1)^{\ell+1}\f$.
 *
 * References:
 * - [1] Erratum: Thermal Casimir effect for Drude metals in the plane-sphere
 *       geometry, Canaguier-Durand, Neto, Lambrecht, Reynaud (2010)
 *       http://journals.aps.org/pra/abstract/10.1103/PhysRevA.83.039905
 * - [2] Thermal Casimir effect for Drude metals in the plane-sphere geometry,
 *       Canaguier-Durand, Neto, Lambrecht, Reynaud (2010),
 *       http://journals.aps.org/pra/abstract/10.1103/PhysRevA.82.012511
 * - [3] Negative Casimir entropies in the plane-sphere geometry, Hartmann, 2014
 * - [4] Casimir effect in the plane-sphere geometry: Beyond the proximity
 *       force approximation, Hartmann, 2018
 *
 * @param [in,out] self CaPS object
 * @param [in] xi_ \f$\xi\mathcal{L}/c\f$
 * @param [in] l angular momentum \f$\ell\f$
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 */
void caps_mie(caps_t *self, double xi_, int l, double *lna, double *lnb)
{
    const double epsilonm1 = caps_epsilonm1_sphere(self, xi_); /* n²-1 */

    if(isinf(epsilonm1))
    {
        /* Mie coefficients for perfect reflectors */
        caps_mie_perf(self, xi_, l, lna, lnb);
        return;
    }

    /* Mie coefficients for arbitrary metals */

    /* χ = ξ*R/(R+L) = ξ/(1+L/R) */
    const double chi    = xi_/(1+self->LbyR);
    const double ln_chi = log(xi_)-log1p(self->LbyR);
    const double ln_l   = logi(l);

    /* Note: n is the refraction index, n_mat the Matsubara index
     * n    = sqrt(ε(ξ,ω_p,γ))
     * ln_n = ln(sqrt(ε)) = ln(ε)/2 = ln(1+(ε-1))/2 = log1p(ε-1)/2
     */
    const double ln_n = log1p(epsilonm1)/2;
    const double n    = exp(ln_n);

    double logIlp, logKlp, logIlm, logKlm, logIlp_nchi, logKlp_nchi, logIlm_nchi, logKlm_nchi;

    bessel_logInKn_half(l,   chi, &logIlp, &logKlp); /* I_{l+0.5}(χ), K_{l+0.5}(χ) */
    bessel_logInKn_half(l-1, chi, &logIlm, &logKlm); /* K_{l-0.5}(χ), K_{l-0.5}(χ) */

    bessel_logInKn_half(l,   n*chi, &logIlp_nchi, &logKlp_nchi); /* I_{l+0.5}(nχ), K_{l+0.5}(nχ) */
    bessel_logInKn_half(l-1, n*chi, &logIlm_nchi, &logKlm_nchi); /* K_{l-0.5}(nχ), K_{l-0.5}(nχ) */

    double ratio   = bessel_ratioI(l-0.5,   chi);
    double ratio_n = bessel_ratioI(l-0.5, n*chi);

    double ln_gammaB = logIlp_nchi + logIlp + ln_chi + log( n*ratio_n-ratio );
    double ln_gammaC = log(epsilonm1)+logIlp_nchi + logadd(ln_chi+logKlm, ln_l+logKlp);
    double ln_gammaD = ln_chi + logadd(logIlp_nchi+logKlm, ln_n+logKlp+logIlm_nchi);

    /* ln_gammaA - ln_gammba_B */
    double num = logIlp_nchi + logIlp + log( -l*epsilonm1 + n*chi*( n*ratio - ratio_n ) );

    *lnb = CAPS_LOGPI-CAPS_LOG2 + ln_gammaB-ln_gammaD;
    *lna = CAPS_LOGPI-CAPS_LOG2 + num - logadd(ln_gammaC, ln_gammaD);
}

/**
 * @brief Calculate Fresnel coefficients \f$r_\mathrm{TE}\f$ and \f$r_\mathrm{TM}\f$ for arbitrary metals
 *
 * This function calculates the Fresnel coefficients \f$r_p = r_p(i\xi, k)\f$
 * for \f$p=\mathrm{TE},\mathrm{TM}\f$.
 *
 * @param [in]     self  CaPS object
 * @param [in]     xi_   \f$\xi\mathcal{L}/c\f$
 * @param [in]     k_    \f$k\mathcal{L}\f$
 * @param [in,out] r_TE  Fresnel coefficient for TE mode
 * @param [in,out] r_TM  Fresnel coefficient for TM mode
 */
void caps_fresnel(caps_t *self, double xi_, double k_, double *r_TE, double *r_TM)
{
    const double epsilonm1 = caps_epsilonm1_plate(self, xi_);

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
     *     β = sqrt( 1 + xi_²/(xi_²+k_²)*(ε-1) ) = sqrt(1+x),
     * where
     *     x = xi_²/(xi_²+k_²)*(ε-1).
     *
     * We calculate x. If x is small, β≈1 and a loss of significance occures
     * when calculating 1-β. For this reason we use sqrtpm1 which calculates
     * β-1.
     */
    const double x = POW_2(xi_)/(POW_2(xi_)+POW_2(k_))*epsilonm1;
    const double betam1 = sqrtpm1(x); /* β-1 */
    *r_TE = -betam1/(2+betam1);
    *r_TM = (epsilonm1-betam1)/(epsilonm1+2+betam1);
}

/*@}*/


/*@}*/


/**
* @name Kernels
*/
/*@{*/

/**
 * @brief Initialize caps_M_t object
 *
 * This object contains all information necessary to compute the matrix
 * elements of the round-trip operator \f$\mathcal{M}^{(m)}(\xi)\f$. It also
 * contains a cache for the Mie coefficients.
 *
 * The returned object can be given to \ref caps_kernel_M to compute the
 * matrix elements of the round-trip operator.
 *
 * @param [in] caps CaPS object
 * @param [in] m azimuthal quantum number \f$m\f$
 * @param [in] xi_ \f$\xi\mathcal{L}/c\f$
 * @retval obj caps_M_t object that can be given to \ref caps_kernel_M
 */
caps_M_t *caps_M_init(caps_t *caps, int m, double xi_)
{
    TERMINATE(xi_ <= 0, "Matsubara frequency must be positive");

    size_t lmin, lmax;
    const int ldim = caps->ldim;
    caps_M_t *self = xmalloc(sizeof(caps_M_t));

    caps_estimate_lminmax(caps, m, &lmin, &lmax);

    self->caps = caps;
    self->m = m;
    self->lmin = lmin;
    self->ldim = ldim;
    self->integration = caps_integrate_init(caps, xi_, m, caps->epsrel);
    self->integration_plasma = NULL;
    self->xi_ = xi_;
    self->al = xmalloc(ldim*sizeof(double));
    self->bl = xmalloc(ldim*sizeof(double));

    for(int j = 0; j < ldim; j++)
        self->al[j] = self->bl[j] = NAN;

    return self;
}

/**
 * @brief Kernel of round-trip matrix
 *
 * This function returns the matrix elements of the round-trip operator
 * \f$\mathcal{M}^{(m)}\f$.
 *
 * The round-trip matrix is a \f$2\ell_\mathrm{dim} \times 2\ell_\mathrm{dim}\f$
 * matrix, the matrix elements start at 0, i.e. \f$0 \le i,j < 2\ell_\mathrm{dim}\f$.
 *
 * This function is intended to be passed as a callback to \ref kernel_logdet. If you
 * want to compute matrix elements of the round-trip operator, it is probably simpler
 * to use \ref caps_M_elem.
 *
 * @param [in] i row
 * @param [in] j column
 * @param [in] args_ caps_M_t object, see \ref caps_M_init
 * @retval Mij \f$\mathcal{M}^{(m)}_{ij}(\xi)\f$
 */
double caps_kernel_M(int i, int j, void *args_)
{
    caps_M_t *args = (caps_M_t *)args_;
    const int lmin = args->lmin;

    #if 1
    /* variant A: (faster)
     * round-trip matrix consists of many 2x2 blocks */
    const int l1 = lmin+i/2, l2 = lmin+j/2;

    char polarizations[2] = { 'E', 'M' };
    char p1 = polarizations[i % 2], p2 = polarizations[j % 2];
    #else
    /* variant B: (slower)
     * round-trip matrix consists of 4 big polarization blocks, EE, EM, ME, MM */
    const int ldim = args->ldim;
    const int l1 = lmin + (i % ldim), l2 = lmin + (j % ldim);

    char p1 = 'E', p2 = 'E';
    if(i >= ldim)
        p1 = 'M';
    if(j >= ldim)
        p2 = 'M';
    #endif

    return caps_M_elem(args, l1, l2, p1, p2);
}

/**
 * @brief Compute matrix elements of round-trip operator
 *
 * This function computes matrix elements of the round-trip operator.
 *
 * Warning: Make sure that lmin <= l1,l2 <= lmax or otherwise the behavior of
 * this function is undefined. You can get lmin and lmax using \ref
 * caps_estimate_lminmax.
 *
 * @param [in] self caps_M_t object, see \ref caps_M_init
 * @param [in] l1 angular momentum \f$\ell_1\f$
 * @param [in] l2 angular momentum \f$\ell_2\f$
 * @param [in] p1 polarization \f$p_1\f$ (E or M)
 * @param [in] p2 polarization \f$p_2\f$ (E or M)
 * @retval elem \f$\mathcal{M}^{(m)}_{\ell_1,\ell_2}(p_1,p_2)\f$
 */
double caps_M_elem(caps_M_t *self, int l1, int l2, char p1, char p2)
{
    const double xi_ = self->xi_; /* xi_ = xi*(L+R)/c */
    const int lmin = self->lmin;
    caps_t *caps = self->caps;
    integration_t *integration = self->integration;

    if(isnan(self->al[l1-lmin]))
        caps_mie(caps, xi_, l1, &self->al[l1-lmin], &self->bl[l1-lmin]);

    if(isnan(self->al[l2-lmin]))
        caps_mie(caps, xi_, l2, &self->al[l2-lmin], &self->bl[l2-lmin]);

    const double lnLambda = caps_lnLambda(l1,l2,self->m);
    const double al1 = self->al[l1-lmin], bl1 = self->bl[l1-lmin];
    const double al2 = self->al[l2-lmin], bl2 = self->bl[l2-lmin];

    if(p1 == p2) /* EE or MM */
    {
        if(p1 == 'E') /* EE */
        {
            sign_t signA_TE, signB_TM;
            double log_A_TE = caps_integrate_A(integration, l1, l2, TE, &signA_TE);
            double log_B_TM = caps_integrate_B(integration, l1, l2, TM, &signB_TM);

            /* √(a_l1*a_l2)*(B_TM - A_TE) */
            double mie = (al1+al2)/2;
            return exp(log_B_TM+mie+lnLambda)*signB_TM-exp(log_A_TE+mie+lnLambda)*signA_TE;
        }
        else /* MM */
        {
            sign_t signA_TM, signB_TE;
            double log_A_TM = caps_integrate_A(integration, l1, l2, TM, &signA_TM);
            double log_B_TE = caps_integrate_B(integration, l1, l2, TE, &signB_TE);

            /* √(b_l1*b_l2)*(A_TM - B_TE) */
            double mie = (bl1+bl2)/2;
            return exp(log_A_TM+mie+lnLambda)*signA_TM-exp(log_B_TE+mie+lnLambda)*signB_TE;
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

            double log_C_TE = caps_integrate_C(integration, l1, l2, TE, &signC_TE);
            double log_D_TM = caps_integrate_D(integration, l1, l2, TM, &signD_TM);

            /* D_TM - C_TE */
            double mie1 = (al1+bl2)/2;
            return exp(log_D_TM+mie1+lnLambda)*signD_TM-exp(log_C_TE+mie1+lnLambda)*signC_TE;
        }
        else /* ME */
        {
            sign_t signC_TM, signD_TE;

            double log_C_TM = caps_integrate_C(integration, l1, l2, TM, &signC_TM);
            double log_D_TE = caps_integrate_D(integration, l1, l2, TE, &signD_TE);

            /* C_TM - D_TE */
            const double mie2 = (bl1+al2)/2;
            return exp(log_C_TM+mie2+lnLambda)*signC_TM-exp(log_D_TE+mie2+lnLambda)*signD_TE;
        }
    }
}

/**
 * @brief Free caps_M_t object
 *
 * Frees memory allocated by \ref caps_M_init.
 *
 * @param [in,out] self caps_M_t object
 */
void caps_M_free(caps_M_t *self)
{
    xfree(self->al);
    xfree(self->bl);
    caps_integrate_free(self->integration);
    xfree(self);
}

/** @brief Kernel for EE block
 *
 * Function that returns matrix elements of the round-trip matrix
 * \f$\mathcal{M}\f$ for \f$\xi=0\f$ and polarization \f$p_1=p_2=\mathrm{E}\f$.
 *
 * See also \ref caps_logdetD0.
 *
 * @param [in] i row (starting from 0)
 * @param [in] j column (starting from 0)
 * @param [in] args_ pointer to caps_M_t object
 * @retval Mij matrix element
 */
double caps_kernel_M0_EE(int i, int j, void *args_)
{
    caps_M_t *args = (caps_M_t *)args_;
    const double y = args->caps->y;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    return exp( (l1+l2+1)*y + lfac(l1+l2) - 0.5*(lfac(l1+m)+lfac(l1-m) + lfac(l2+m)+lfac(l2-m)) );
}

/** @brief Kernel for MM block (plasma model)
 *
 * Function that returns matrix elements of round-trip matrix \f$\mathcal{M}\f$
 * for \f$\xi=0\f$ and polarization \f$p_1=p_2=\mathrm{M}\f$ (plasma model).
 *
 * See also \ref caps_logdetD0.
 *
 * @param [in] i row (starting from 0)
 * @param [in] j column (starting from 0)
 * @param [in] args_ pointer to caps_M_t object
 * @retval Mij matrix element
 */
double caps_kernel_M0_MM_plasma(int i, int j, void *args_)
{
    caps_M_t *args = (caps_M_t *)args_;

    /* geometry */
    const double y = args->caps->y;

    integration_plasma_t *integration_plasma = args->integration_plasma;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    /* I: value of integral
     * ratio1 = I_{l1-1/2}/I_{l1+1/2}
     * ratio2 = I_{l2-1/2}/I_{l2+1/2}
     */
    double ratio1,ratio2;
    double I = caps_integrate_plasma(integration_plasma, l1, l2, m, &ratio1, &ratio2);

    double alpha = integration_plasma->alpha;
    double factor1 = 1-(2*l1+1.)/(alpha*ratio1);
    double factor2 = 1-(2*l2+1.)/(alpha*ratio2);
    double factor = sqrt(l1*factor1/(l1+1.) * l2*factor2/(l2+1.)) ;

    return factor*exp( lfac(l1+l2)-0.5*(lfac(l1+m)+lfac(l1-m)+lfac(l2+m)+lfac(l2-m)) + (l1+l2+1)*y )*I;
}

/** @brief Kernel for MM block
 *
 * Function that returns matrix elements of round-trip matrix \f$\mathcal{M}\f$
 * for \f$\xi=0\f$ and polarization \f$p_1=p_2=\mathrm{M}\f$.
 *
 * See also \ref caps_logdetD0.
 *
 * @param [in] i row (starting from 0)
 * @param [in] j column (starting from 0)
 * @param [in] args_ pointer to caps_M_t object
 * @retval Mij matrix element
 */
double caps_kernel_M0_MM(int i, int j, void *args_)
{
    caps_M_t *args = (caps_M_t *)args_;
    const double y = args->caps->y;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    return exp( (l1+l2+1)*y + lfac(l1+l2) - 0.5*(lfac(l1+m)+lfac(l1-m) + lfac(l2+m)+lfac(l2-m)) )*sqrt(((double)l1*(double)l2)/((l1+1.)*(l2+1.)));
}

/*@}*/

/**
 * @name Compute determinants
 */
/*@{*/

/** @brief Compute \f$\log\det\mathcal{D}^{(m)}\left(\frac{\xi\mathcal{L}}{c}\right)\f$
 *
 * This function computes the logarithm of the determinant of the scattering
 * matrix for the frequency \f$\xi\mathcal{L}/c\f$ and quantum number
 * \f$m\f$.
 *
 * For \f$\xi=0\f$ see \ref caps_logdetD0.
 *
 * @param self CaPS object
 * @param xi_ \f$\xi\mathcal{L}/c > 0\f$
 * @param m quantum number \f$m\f$
 * @retval logdetD
 */
double caps_logdetD(caps_t *self, double xi_, int m)
{
    TERMINATE(xi_ <= 0, "Matsubara frequency must be positive");

    const int sym_spd = 2; /* matrix is symmetric and positive definite */
    const int dim = 2*self->ldim;

    caps_M_t *args = caps_M_init(self, m, xi_);
    double logdet = kernel_logdet(dim, &caps_kernel_M, args, sym_spd, self->detalg);
    caps_M_free(args);

    return logdet;
}


/** @brief Compute \f$\log\det\mathcal{D}^{(m)}(\xi=0)\f$ for EE and/or MM contribution
 *
 * Compute numerically for a given value of \f$m\f$ the contribution of the
 * polarization block EE and/or MM. If EE, MM or MM_plasma is NULL, the value
 * will not be computed.
 *
 * For Drude metals there exists an analytical formula to compute logdetD, see
 * \ref caps_ht_drude.
 *
 * For perfect reflectors see also \ref caps_ht_perf.
 *
 * For the Plasma model see also \ref caps_ht_plasma.
 *
 * @param [in]  self CaPS object
 * @param [in]  m quantum number \f$m\f$
 * @param [in]  omegap plasma frequency in rad/s (only used to compute MM_plasma)
 * @param [out] EE pointer to store contribution for EE block
 * @param [out] MM pointer to store contribution for MM block
 * @param [out] MM_plasma pointer to store contribution for MM block (Plasma model)
 */
void caps_logdetD0(caps_t *self, int m, double omegap, double *EE, double *MM, double *MM_plasma)
{
    size_t lmin, lmax;
    const int sym_spd = 2; /* matrix is symmetric and positive definite */
    const int ldim = self->ldim;
    detalg_t detalg = self->detalg;

    caps_estimate_lminmax(self, m, &lmin, &lmax);

    caps_M_t args = {
        .caps = self,
        .m = m,
        .xi_ = 0,
        .integration = NULL,
        .integration_plasma = NULL,
        .al = NULL,
        .bl = NULL,
        .lmin = lmin
    };

    if(EE != NULL)
        *EE = kernel_logdet(ldim, &caps_kernel_M0_EE, &args, sym_spd, detalg);

    if(MM != NULL)
        *MM = kernel_logdet(ldim, &caps_kernel_M0_MM, &args, sym_spd, detalg);

    if(MM_plasma != NULL)
    {
        args.integration_plasma = caps_integrate_plasma_init(self, omegap, self->epsrel);
        *MM_plasma = kernel_logdet(ldim, &caps_kernel_M0_MM_plasma, &args, sym_spd, detalg);
        caps_integrate_plasma_free(args.integration_plasma);
    }
}

/*@}*/

/**
* @name high-temperature limit
*/
/*@{*/

/** @brief Compute high-temperature limit for Drude metals
 *
 * For Drude metals the Fresnel coefficients become \f$r_\mathrm{TM}=1\f$,
 * \f$r_\mathrm{TE}=0\f$ for \f$\xi\to 0\f$, i.e. only the EE polarization
 * block needs to be considered.
 *
 * For Drude the free energy for \f$\xi=0\f$ can be computed analytically. We
 * use Eq. (8) from Ref. [1] to compute the contribution.
 *
 * References:
 *  - [1] Bimonte, Emig, "Exact results for classical Casimir interactions:
 *    Dirichlet and Drude model in the sphere-sphere and sphere-plane geometry",
 *    Phys. Rev. Lett. 109 (2012), https://doi.org/10.1103/PhysRevLett.109.160403
 *
 * @param [in] caps CaPS object
 * @retval F free energy in units of \f$k_\mathrm{B}T\f$
 */
double caps_ht_drude(caps_t *caps)
{
    const double x  = caps->LbyR; /* L/R */
    const double x2 = POW_2(x);      /* (L/R)² */
    const double Z  = (1+x)*(1-sqrt((x2+2*x)/(x2+2*x+1)));

    double sum1 = 0, sum2 = 0;
    for(int l = 1; true; l++)
    {
        const double a = pow(Z, 2*l); /* Z^(2l) */
        const double b = a*Z;         /* Z^(2l+1) */
        const double v1 = (2*l+1)*log1p(-b);
        const double v2 = b*(1-a)/(1-b);

        sum1 += v1;
        sum2 += v2;

        if(fabs(v1/sum1) < 1e-15 && fabs(v2/sum2) < 1e-15)
            break;
    }

    return 0.5*(sum1 + log1p(-(1-POW_2(Z))*sum2));
}

/** @brief Compute free energy in the high-temperature limit for perfect reflectors
 *
 * For perfect reflectors the Fresnel coefficients become
 * \f$r_\mathrm{TM}=1\f$, \f$r_\mathrm{TE}=-1\f$ in the limit \f$\xi\to 0\f$,
 * and only the polarization blocks EE and MM need to be considered.
 *
 * The contribution for EE, i.e. Drude, can be computed analytically, see \ref
 * caps_ht_drude. For the MM block we numerically compute the determinants
 * up to \f$m = M\f$ until
 * \f[
 *      \frac{\log\det\mathcal{D}^{(M)}(0)}{{\sum_{m=0}^M}^\prime \log\det\mathcal{D}^{(m)}(0)} < \epsilon \,.
 * \f]
 * We use Ref. [1] to compute the contribution for \f$m = 0\f$.
 *
 * References:
 *  - [1] Bimonte, Classical Casimir interaction of perfectly conducting sphere
 *    and plate (2017), https://arxiv.org/abs/1701.06461
 *
 * @param [in] caps CaPS object
 * @param [in] eps \f$\epsilon\f$ abort criterion
 * @retval energy free energy in units of \f$k_\mathrm{B}T\f$
 */
double caps_ht_perf(caps_t *caps, double eps)
{
    const double LbyR = caps->LbyR;
    const double drude = caps_ht_drude(caps);

    double MM = 0;

    /* m = 0, see Ref. [1] */
    const double Z = 1./(1+LbyR+sqrt(LbyR*(2+LbyR))); /* Eq. (7) */
    for(int l = 0; true; l++)
    {
        const double v = 0.5*log1p(-pow(Z,2*l+3)); /* Eq. (17) */
        MM += v;
        if(fabs(v/MM) < eps)
            break;
    }

    for(int m = 1; true; m++)
    {
        double v;
        caps_logdetD0(caps, m, 0, NULL, &v, NULL);

        MM += v;

        if(fabs(v/MM) < eps)
            return drude+MM;
    }
}

/** @brief Compute free energy in the high-temperature limit for plasma model
 *
 * The abort criterion eps is the same as in \ref caps_ht_perf.
 *
 * @param [in] caps CaPS object
 * @param [in] omegap plasma frequency in rad/s
 * @param [in] eps abort criterion
 * @retval F free energy in units of \f$k_\mathrm{B}T\f$
 */
double caps_ht_plasma(caps_t *caps, double omegap, double eps)
{
    const double drude = caps_ht_drude(caps);
    double MM_plasma = 0;

    for(int m = 0; true; m++)
    {
        double v;
        caps_logdetD0(caps, m, omegap, NULL, NULL, &v);

        if(m == 0)
            MM_plasma += v/2;
        else
            MM_plasma += v;

        if(fabs(v/MM_plasma) < eps)
            return drude+MM_plasma;
    }
}

/*@}*/
