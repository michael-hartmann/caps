/**
 * @file   libcasimir.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   August, 2017
 * @brief  library to calculate the free Casimir energy in the plane-sphere geometry
 */

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "quadpack.h"
#include "constants.h"
#include "bessel.h"
#include "integration.h"
#include "libcasimir.h"
#include "matrix.h"
#include "logfac.h"
#include "misc.h"
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
 * @param [in]  size length of str
 * @retval retval bytes written if successful
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

    fprintf(stream, "%sL/R    = %.16g\n", prefix, self->LbyR);
    fprintf(stream, "%sldim   = %d\n",    prefix, self->ldim);
    fprintf(stream, "%sepsrel = %g\n",    prefix, self->epsrel);
    fprintf(stream, "%sdetalg = %s\n",    prefix, detalg_str);
}
/*@}*/

/**
* @name various functions
*/
/*@{*/

/**
 * @brief Calculate logarithm \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$
 *
 * This function returns the logarithm of \f$\Lambda_{\ell_1 \ell_2}^{(m)}\f$ for
 * \f$\ell_1,\ell_2,m\f$.
 * \f[
 *      \Lambda_{\ell_1,\ell_2}^{(m)} = \frac{2 N_{\ell_1,m} N_{\ell_2,m}}{\sqrt{\ell_1 (\ell_1+1) \ell_2 (\ell_2+1)}}
 * \f]
 *
 * The values are computed using the lgamma function to avoid overflows.
 *
 * Symmetries: \f$\Lambda_{\ell_1,\ell_2}^{(m)} = \Lambda_{\ell_2,\ell_1}^{(m)}\f$
 *
 * @param [in]  l1 l1>0
 * @param [in]  l2 l2>0
 * @param [in]  m  m <= l1 and m <= l2
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
 * This function calculates the Fresnel coefficients for TE and TM mode.
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
     * We calculate x. If x is small, β≈1 and a loss of significance occures
     * when calculating 1-β. For this reason we use sqrtpm1 which calculates
     * β-1.
     *
     * Note: ξ=nT
     */
    const double x = pow_2(nT)/(pow_2(nT)+pow_2(k))*epsilonm1;
    const double betam1 = sqrtpm1(x); /* β-1 */
    *r_TE = -betam1/(2+betam1);
    *r_TM = (epsilonm1-betam1)/(epsilonm1+2+betam1);
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
    self->y = -M_LOG2-log1p(LbyR); /* log( (R/(L+R))/2 ) */

    /* dimension of vector space */
    self->ldim = ceil(MAX(CASIMIR_MINIMUM_LDIM, CASIMIR_FACTOR_LDIM/LbyR));

    /* perfect reflectors */
    self->epsilonm1 = casimir_epsilonm1_perf;
    self->rp        = casimir_rp;
    self->lnab      = casimir_lnab;
    self->userdata  = NULL;

    /* relative error for integration */
    self->epsrel = CASIMIR_EPSREL;

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
    xfree(self);
}

/**
 * @brief Set relative rror for numerical integration
 *
 * Set relative error for numerical integration.
 *
 * @param [in] self Casimir object
 * @param [in] epsrel relative error
 * @retval 0 if an error occured
 * @retval 1 on success
 */
int casimir_set_epsrel(casimir_t *self, double epsrel)
{
    if(epsrel <= 0)
        return 0;

    self->epsrel = epsrel;
    return 1;
}

/**
 * @brief Get relative error for numerical integration
 *
 * @retval epsrel relative error
 */
double casimir_get_epsrel(casimir_t *self)
{
    return self->epsrel;
}

/**
 * @brief Set material parameters for generic metals
 *
 * Set dielectric function of material.
 *
 * The Fresnel coefficient and the Mie coefficient depend on the dielectric
 * function epsilon(i*xi). By default, perfect reflectors with a dielectric
 * function epsilon(i*xi)=inf are used.
 *
 * However, you can also specify an arbitrary function for epsilon(i*xi).
 * userdata is an arbitrary pointer that will be given to the callback
 * function.
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

void casimir_set_rp(casimir_t *self, void (*rp)(struct casimir *self, double nT, double k, double *r_TE, double *r_TM))
{
    self->rp = rp;
}

void casimir_set_lnab(casimir_t *self, void (*lnab)(struct casimir *self, double nT, int l, double *lna, double *lnb))
{
    self->lnab = lnab;
}

/**
 * @brief Set algorithm to calculate deterimant
 *
 * The algorithm is given by detalg. Make sure that detalg contains a valid
 * algorithm, otherwise the computation will print a warning on runtime and
 * default to DETALG_LU.
 *
 * detalg may be: DETALG_LU, DETALG_HODLR
 *
 * @param [in,out] self Casimir object
 * @param [in] detalg algorithm to compute determinant
 */
void casimir_set_detalg(casimir_t *self, detalg_t detalg)
{
    if(detalg != DETALG_LU)
        self->detalg = detalg;
#ifdef SUPPORT_LAPACK
    else
        self->detalg = detalg;
#endif
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

/*@}*/


/**
 * @name Mie coefficients
 */
/*@{*/

/**
 * @brief Calculate Mie coefficients for perfect reflectors
 *
 * This function calculates the logarithms of the Mie coefficients
 * \f$a_\ell(i\chi)\f$ and \f$b_\ell(i\chi)\f$ for perfect reflectors. The Mie
 * coefficients are evaluated at the argument \f$\chi=nT R/(R+L)\f$.
 *
 * The signs are given by \f$\mathrm{sgn}(a_\ell) = (-1)^\ell\f$, \f$\mathrm{sgn}(b_\ell) = (-1)^{\ell+1}\f$.
 *
 * lna and lnb must be valid pointers and must not be NULL.
 *
 * Restrictions: \f$\ell \ge 1\f$, \f$\ell \ge 0\f$
 *
 * @param [in,out] self Casimir object
 * @param [in] nT Matsubara frequency
 * @param [in] l angular momentum \f$\ell\f$
 * @param [out] ln_a logarithm of \f$a_\ell\f$
 * @param [out] ln_b logarithm of \f$b_\ell\f$
 */
void casimir_lnab_perf(casimir_t *self, double nT, int l, double *lna, double *lnb)
{
    double lnKlp,lnKlm,lnIlm,lnIlp;
    const double chi = nT/(1+self->LbyR); /* xi*R/(R+L) = xi/(1+L/R) */
    const double log_chi = log(nT)-log1p(self->LbyR);

    /* we could do both calculations together. but it doesn't cost much time -
     * so why bother?
     */
    bessel_lnInuKnu(l-1, chi, &lnIlm, &lnKlm);
    bessel_lnInuKnu(l,   chi, &lnIlp, &lnKlp);

    /* Calculate b_l(chi), i.e. lnb */
    *lnb = M_LOGPI-M_LOG2+lnIlp-lnKlp;

    /* We want to calculate
     *
     * b_l(χ) = (-1)^(l+1)* π/2 * Ip/Im
     * a_l(χ) = (-1)^l    * π/2 * ( χ*Ilm - l*Ilp )/( l*Kp + χ*Km )
     *        = (-1)^l    * π/2 * Ilp * ( χ*Ilm/Ilp - l )/( l*Kp + χ*Km )
     *                                      \-----/
     *                                       ratio
     *
     * where Ip = I_{l+1/2}(χ), Im = I_{l-1/2}(χ), and similar for Kp and Km.
     *
     * Also note that all terms in brackets are positive.
     */

    /* numerator and denominator to calculate al */
    double ratio = bessel_continued_fraction(l-1, chi);
    double numerator = M_LOGPI-M_LOG2 + lnIlp + log(chi*ratio - l);
    double denominator = logadd(logi(l)+lnKlp, log_chi+lnKlm);

    *lna = numerator-denominator;
}

/**
 * @brief Return logarithm of Mie coefficients \f$a_\ell\f$, \f$b_\ell\f$ for arbitrary metals
 *
 * For \f$\omega_\mathrm{P} = \infty\f$ the Mie coefficient for perfect
 * reflectors are returned (see \ref casimir_lnab_perf).
 *
 * lna, lnb, sign_a and sign_b must be valid pointers and must not be NULL.
 *
 * For generic metals we calculate the Mie coefficients \f$a_\ell(i\xi)\f$ und \f$b_\ell(i\xi)\f$ using
 * the expressions taken from [1]. Ref. [1] is the erratum to [2]. Please note
 * that the equations (3.30) and (3.31) in [3] are wrong.
 *
 * Note: If sla =~ slb or slc =~ sld, there is a loss of significance when
 * calculating sla-slb or slc-sld.
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
 *
 * @param [in,out] self Casimir object
 * @param [in] nT Matsubara frequency
 * @param [in] l angular momentum \f$\ell\f$
 * @param [out] lna logarithm of Mie coefficient \f$a_\ell\f$
 * @param [out] lnb logarithm of Mie coefficient \f$b_\ell\f$
 */
void casimir_lnab(casimir_t *self, double nT, int l, double *lna, double *lnb)
{
    /* ξ = nT */
    const double epsilonm1 = casimir_epsilonm1(self, nT); /* n²-1 */

    if(isinf(epsilonm1))
    {
        /* Mie coefficients for perfect reflectors */
        casimir_lnab_perf(self, nT, l, lna, lnb);
        return;
    }

    /* Mie coefficients for arbitrary metals */

    /* χ = ξ*R/(R+L) = ξ/(1+L/R) */
    const double chi    = nT/(1+self->LbyR);
    const double ln_chi = log(nT)-log1p(self->LbyR);
    const double ln_l   = logi(l);

    /* Note: n is the refraction index, n_mat the Matsubara index
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

    double ratio   = bessel_continued_fraction(l-1,   chi);
    double ratio_n = bessel_continued_fraction(l-1, n*chi);

    double ln_gammaB = lnIlp_nchi + lnIlp + ln_chi + log( n*ratio_n-ratio );
    double ln_gammaC = log(epsilonm1)+lnIlp_nchi + logadd(ln_chi+lnKlm, ln_l+lnKlp);
    double ln_gammaD = ln_chi + logadd(lnIlp_nchi+lnKlm, ln_n+lnKlp+lnIlm_nchi);

    /* ln_gammaA - ln_gammba_B */
    double num = lnIlp_nchi + lnIlp + log( -l*epsilonm1 + n*chi*( n*ratio - ratio_n ) );

    *lnb = M_LOGPI-M_LOG2 + ln_gammaB-ln_gammaD;
    *lna = M_LOGPI-M_LOG2 + num - logadd(ln_gammaC, ln_gammaD);
}
/*@}*/


/* XXX simplify XXX */
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
    TERMINATE(nT <= 0, "Matsubara frequency must be positive");

    size_t lmin, lmax;
    const int ldim = casimir->ldim;
    casimir_M_t *self = xmalloc(sizeof(casimir_M_t));

    casimir_estimate_lminmax(casimir, m, &lmin, &lmax);

    self->casimir = casimir;
    self->m = m;
    self->lmin = lmin;
    self->integration = casimir_integrate_init(casimir, nT, m, casimir->epsrel);
    self->integration_plasma = NULL;
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
        casimir->lnab(casimir, nT, l1, &self->al[l1-lmin], &self->bl[l1-lmin]);

    if(isnan(self->al[l2-lmin]))
        casimir->lnab(casimir, nT, l2, &self->al[l2-lmin], &self->bl[l2-lmin]);

    const double al1 = self->al[l1-lmin], bl1 = self->bl[l1-lmin];
    const double al2 = self->al[l2-lmin], bl2 = self->bl[l2-lmin];

    if(p1 == p2) /* EE or MM */
    {
        if(p1 == 'E') /* EE */
        {
            sign_t signA_TE, signB_TM;
            double log_A_TE = casimir_integrate_A(integration, l1, l2, TE, &signA_TE);
            double log_B_TM = casimir_integrate_B(integration, l1, l2, TM, &signB_TM);

            /* √(a_l1*a_l2)*(B_TM - A_TE) */
            double mie = (al1+al2)/2;
            return exp(log_B_TM+mie)*signB_TM-exp(log_A_TE+mie)*signA_TE;
        }
        else /* MM */
        {
            sign_t signA_TM, signB_TE;
            double log_A_TM = casimir_integrate_A(integration, l1, l2, TM, &signA_TM);
            double log_B_TE = casimir_integrate_B(integration, l1, l2, TE, &signB_TE);

            /* √(b_l1*b_l2)*(A_TM - B_TE) */
            double mie = (bl1+bl2)/2;
            return exp(log_A_TM+mie)*signA_TM-exp(log_B_TE+mie)*signB_TE;
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

            /* D_TM - C_TE */
            double mie1 = (al1+bl2)/2;
            return exp(log_D_TM+mie1)*signD_TM-exp(log_C_TE+mie1)*signC_TE;
        }
        else /* ME */
        {
            sign_t signC_TM, signD_TE;

            double log_C_TM = casimir_integrate_C(integration, l1, l2, TM, &signC_TM);
            double log_D_TE = casimir_integrate_D(integration, l1, l2, TE, &signD_TE);

            /* C_TM - D_TE */
            const double mie2 = (bl1+al2)/2;
            return exp(log_C_TM+mie2)*signC_TM-exp(log_D_TE+mie2)*signD_TE;
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


/** @brief Compute log det D^m(xi)
 *
 * This function will compute the logarithm of the determinant of the
 * scattering matrix for Matsubara frequency nT and quantum number m.
 *
 * Either LU decomposition (slow) or method for HODLR matrices (fast) will be
 * used, see \ref casimir_set_detalg.
 *
 * Please note that the scaled Matsubara frequency must be positive.
 *
 * @param self Casimir object
 * @param nT Matsubara frequency
 * @param m quantum number m
 * @retval logdetD(xi,m)
 */
double casimir_logdetD(casimir_t *self, double nT, int m)
{
    TERMINATE(nT <= 0, "Matsubara frequency must be positive");

    const int is_symmetric = 1;
    const int dim = 2*self->ldim;

    casimir_M_t *args = casimir_M_init(self, m, nT);
    double logdet = kernel_logdet(dim, &casimir_kernel_M, args, is_symmetric, self->detalg);
    casimir_M_free(args);

    return logdet;
}

/*@}*/


/**
 * @name Matsubara frequency nT=0
 */
/*@{*/

/** @brief Compute logdet D(xi=0) for Drude metals
 *
 * For Drude metals the Fresnel coefficients become r_TM=1, r_TE=0 for xi->0,
 * i.e. only the EE polarization block needs to be considered.
 *
 * For Drude the free energy for the xi=0 can be computed analytically. We use
 * Eq. (8) from Ref. [1] to compute the contribution.
 *
 * References:
 *  - [1] Bimonte, Emig, "Exact results for classical Casimir interactions:
 *    Dirichlet and Drude model in the sphere-sphere and sphere-plane geometry",
 *    Phys. Rev. Lett. 109 (2012), https://doi.org/10.1103/PhysRevLett.109.160403
 *
 * @param [in] casimir Casimir object
 * @retval logdetD log det D(xi=0) for Drude metals
 */
double casimir_logdetD0_drude(casimir_t *casimir)
{
    const double x  = casimir->LbyR; /* L/R */
    const double x2 = pow_2(x);      /* (L/R)² */
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

    return 0.5*(sum1 + log1p(-(1-pow_2(Z))*sum2));
}


/** @brief Compute logdet D(xi=0) for perect conductors
 *
 * For Drude metals the Fresnel coefficients are r_TM=1, r_TE=-1. In the limit
 * xi->0 only the polarization blocks EE and MM need to be considered.
 *
 * The contribution for EE, i.e. Drude, can be computed analytically, see \ref
 * casimir_logdetD0_drude. For the MM block we numerically compute the
 * determinants of the m > 0 contributions until
 *      logdetD^m/logdetD^(m=0) < eps.
 * We use Ref. [1] to compute the contribution for m = 0.
 *
 * References:
 *  - [1] Bimonte, Classical Casimir interaction of perfectly conducting sphere
 *    and plate (2017), https://arxiv.org/abs/1701.06461
 *
 * @param [in] casimir Casimir object
 * @param [in] eps accuracy
 * @retval logdetD log det D(xi=0) for perfect conductors
 */
double casimir_logdetD0_pc(casimir_t *casimir, double eps)
{
    const double LbyR = casimir->LbyR;
    const double drude = casimir_logdetD0_drude(casimir);

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
        casimir_logdetD0(casimir, m, 0, NULL, &v, NULL);

        MM += v;

        if(fabs(v/MM) < eps)
            return drude+MM;
    }
}


/** @brief Compute logdetD for xi=0 for EE and/or MM contribution
 *
 * Compute numerically for a given value of m the contribution of the
 * polarization block EE and/or MM. If logdet_EE or logdet_MM is NULL, the
 * value is not computed.
 *
 * The EE block corresponds to Drude metals. For Drude metals there exists an
 * analytical formula to compute logdetD, see \ref casimir_logdetD0_drude.
 *
 * @param [in]  self Casimir object
 * @param [in]  m quantum number m
 * @param [out] EE pointer to store contribution for EE block
 * @param [out] MM pointer to store contribution for MM block
 */
void casimir_logdetD0(casimir_t *self, int m, double omegap, double *EE, double *MM, double *MM_plasma)
{
    size_t lmin, lmax;
    const int is_symmetric = 1, ldim = self->ldim;
    detalg_t detalg = self->detalg;

    casimir_estimate_lminmax(self, m, &lmin, &lmax);

    casimir_M_t args = {
        .casimir = self,
        .m = m,
        .nT = 0,
        .integration = NULL,
        .integration_plasma = NULL,
        .al = NULL,
        .bl = NULL,
        .lmin = lmin
    };

    if(EE != NULL)
        *EE = kernel_logdet(ldim, &casimir_kernel_M0_EE, &args, is_symmetric, detalg);

    if(MM != NULL)
        *MM = kernel_logdet(ldim, &casimir_kernel_M0_MM, &args, is_symmetric, detalg);

    if(MM_plasma != NULL)
    {
        args.integration_plasma = casimir_integrate_plasma_init(self, omegap, self->epsrel);
        *MM_plasma = kernel_logdet(ldim, &casimir_kernel_M0_MM_plasma, &args, is_symmetric, detalg);
        casimir_integrate_plasma_free(args.integration_plasma);
    }
}

/** @brief Kernel for EE block
 *
 * Function that returns matrix elements of round-trip matrix M for xi=0 and
 * polarization p1=p2=E.
 *
 * @param [in] i row (starting from 0)
 * @param [in] j column (starting from 0)
 * @param [in] args pointer to casimir_M_t object
 */
double casimir_kernel_M0_EE(int i, int j, void *args_)
{
    casimir_M_t *args = (casimir_M_t *)args_;
    const double y = args->casimir->y;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    return exp( (l1+l2+1)*y + lfac(l1+l2) - 0.5*(lfac(l1+m)+lfac(l1-m) + lfac(l2+m)+lfac(l2-m)) );
}

double casimir_kernel_M0_MM_plasma(int i, int j, void *args_)
{
    casimir_M_t *args = (casimir_M_t *)args_;

    /* geometry */
    const double y = args->casimir->y;

    integration_plasma_t *integration_plasma = args->integration_plasma;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    /* value of integral */
    /* ratio1 = I_{l1-1/2}/I_{l1+1/2} ; ratio2 = I_{l2-1/2}/I_{l2+1/2} */
    double ratio1,ratio2;
    double I = casimir_integrate_plasma(integration_plasma, l1, l2, m, &ratio1, &ratio2);

    double alpha = integration_plasma->alpha;
    double factor1 = 1-(2*l1+1.)/(alpha*ratio1);
    double factor2 = 1-(2*l2+1.)/(alpha*ratio2);
    double factor = sqrt(l1*factor1/(l1+1.) * l2*factor2/(l2+1.)) ;

    return factor*exp( lfac(l1+l2)-0.5*(lfac(l1+m)+lfac(l1-m)+lfac(l2+m)+lfac(l2-m)) + (l1+l2+1)*y )*I;
}

/** @brief Kernel for MM block
 *
 * Function that returns matrix elements of round-trip matrix M for xi=0 and
 * polarization p1=p2=M.
 *
 * @param [in] i row (starting from 0)
 * @param [in] j column (starting from 0)
 * @param [in] args pointer to casimir_M_t object
 */
double casimir_kernel_M0_MM(int i, int j, void *args_)
{
    casimir_M_t *args = (casimir_M_t *)args_;
    const double y = args->casimir->y;
    const int lmin = args->lmin;
    const int l1 = i+lmin, l2 = j+lmin, m = args->m;

    return exp( (l1+l2+1)*y + lfac(l1+l2) - 0.5*(lfac(l1+m)+lfac(l1-m) + lfac(l2+m)+lfac(l2-m)) )*sqrt(((double)l1*(double)l2)/((l1+1.)*(l2+1.)));
}


double casimir_logdetD0_plasma(casimir_t *casimir, double omegap, double eps)
{
    const double drude = casimir_logdetD0_drude(casimir);
    double MM_plasma = 0;

    for(int m = 0; true; m++)
    {
        double v;
        casimir_logdetD0(casimir, m, omegap, NULL, NULL, &v);

        if(m == 0)
            MM_plasma += v/2;
        else
            MM_plasma += v;

        if(fabs(v/MM_plasma) < eps)
            return drude+MM_plasma;
    }
}

/*@}*/
