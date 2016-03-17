/**
 * @file   integration_drude.c
 * @author Erik Buchenau <e.buchenau@live.de>
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   January, 2016
 * @brief  Perform integration for Drude planes
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "floattypes.h"
#include "utils.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_drude.h"
#include "gausslaguerre.h"
#include "plm_cache.h"

/* **************************************** */
/*
 * Choose one method of integration
 * Romberg-Integration is recommended
 */
//#define INTEGRATION_GAUSS_LAGUERRE
#define INTEGRATION_ROMBERG

/* **************************************** */

/*
 * Tunable parameter for the integration. A small value will gather the points
 * around the maximum of the integrand, while a large value will spread the
 * points.
 */
#define ALPHA 2

/*
 * The accuracy of the Drude-Integration
 */
#define ACCURACY DRUDE_INTEG_ACCURACY


/*
 * Cache the values of the legendre polynomials. If not defined, the program
 * will caclulate the values of plm with the function "plm_PlmPlm" from sfunc.c.
 * ----- Experimental code ------
 */
#define USE_PLM_CACHE



#if (defined INTEGRATION_ROMBERG && defined INTEGRATION_GAUSS_LAGUERRE)
#error Choose only ONE method of Integration
#endif

#if !((defined INTEGRATION_GAUSS_LAGUERRE) || (defined INTEGRATION_ROMBERG))
#error Choose AT LEAST one method of Integration
#endif


struct drude_error {
    float80        lnA, lnB, lnC, lnD;
};


/*
 * If we choose ALPHA = 2, we can use faster functions to calculate x^ALPHA, x^(ALPHA-1)
 * and x^(ALPHA+1). This is important for integrands_drude_u (which is called very often).
 */
#if ALPHA == 2

#define pow80_alpha(x) pow_2(x)
#define pow80_alpha_p(x) pow_3(x)
#define pow80_alpha_m(x) (x)

#else

#define pow80_alpha(x) pow80(x, ALPHA)
#define pow80_alpha_p(x) pow80(x, ALPHA + 1)
#define pow80_alpha_m(x) pow80(x, ALPHA - 1)

#endif


/* ******************** Prototypes ******************** */
#ifdef INTEGRATION_ROMBERG
static inline void integrate_romberg(struct integ_context* context,
                                     integrands_drude_t* result,
                                     struct drude_error* error);
#endif

void integrate_gauss_laguerre(casimir_t *self,
                              casimir_integrals_t *cint,
                              int l1, int l2, int m, int n, double T);

void integrands_drude(float80 x, integrands_drude_t *integrands,
                      struct integ_context* context,
                      unsigned int index,
                      unsigned int iteration);

/* ****************** End Prototypes ****************** */

#ifdef INTEGRATION_ROMBERG
/** @brief Evaluate the Integrands A,B,C and D without any prefactors
 *
 * This function calculates
 *    \frac{r_p e^{-x}}{x^2 + 2\tau x} P_{l1}^m \left(1+\frac{x}{\tau}\right) P_{l1}^m \left(1+\frac{x}{\tau}\right),
 *    r_p e^{-x} (x^2 + 2\tau x) P'_{l1}^m \left(1+\frac{x}{\tau}\right) P'_{l1}^m \left(1+\frac{x}{\tau}\right),
 *    r_p e^{-x} P_{l1}^m \left(1+\frac{x}{\tau}\right) P'_{l2}^m \left(1+\frac{x}{\tau}\right).
 *
 * @param [in]  x          Argument x of the Integrand
 * @param [out] integrands logarithms of values and signs of the integrand
 * @param [in]  self       Casimir object
 * @param [in]  nT         Scaled Temperature
 * @param [in]  l1         Value of l1
 * @param [in]  l2         Value of l2
 * @param [in]  m          Value of m
 */
void integrands_drude(float80 x, integrands_drude_t *integrands,
                      struct integ_context* context,
                      unsigned int index,
                      unsigned int iteration)
{
    plm_combination_t comb;
    const float80 tau = 2 * context->nT;
    const float80 k   = sqrt80(pow_2(x) / 4 + context->nT * x);
    float80 r_TE, r_TM;

    casimir_rp(context->casimir, context->nT, k, &r_TE, &r_TM);
    const float80 lnr_TE = log80(-r_TE);
    const float80 lnr_TM = log80(r_TM);

    if( fabs80(x/tau - 1.0) >= LDBL_EPSILON )
    {
#ifdef USE_PLM_CACHE
        plm_cache_PlmPlm(context, 1.0+x/tau, &comb, index, iteration);
#else
        plm_PlmPlm(context->l1, context->l2, context->m, 1.0+x/tau, &comb);
#endif
    }
    else
    {
        comb.lnPl1mPl2m      = -INFINITY;
        comb.lnPl1mdPl2m     = -INFINITY;
        comb.lndPl1mPl2m     = -INFINITY;
        comb.lndPl1mdPl2m    = -INFINITY;
        comb.sign_Pl1mPl2m   = 1;
        comb.sign_Pl1mdPl2m  = 1;
        comb.sign_dPl1mPl2m  = 1;
        comb.sign_dPl1mdPl2m = 1;
    }
    const float80 log_factor = log80(pow_2(x)+2*tau*x);

    const float80 A = comb.lnPl1mPl2m - log_factor - x;
    integrands->lnA_TE = lnr_TE + A;
    integrands->lnA_TM = lnr_TM + A;
    integrands->sign_A = comb.sign_Pl1mPl2m;

    const float80 B = comb.lndPl1mdPl2m + log_factor - x;
    integrands->lnB_TE = lnr_TE + B;
    integrands->lnB_TM = lnr_TM + B;
    integrands->sign_B = comb.sign_dPl1mdPl2m;

    const float80 C = comb.lnPl1mdPl2m - x;
    integrands->lnC_TE = lnr_TE + C;
    integrands->lnC_TM = lnr_TM + C;
    integrands->sign_C = comb.sign_Pl1mdPl2m;

    const float80 D = comb.lndPl1mPl2m;
    integrands->lnD_TE = lnr_TE + D - x;
    integrands->lnD_TM = lnr_TM + D - x;
    integrands->sign_D = comb.sign_dPl1mPl2m;
}


/** @brief Evaluate the Integrands A,B,C,D without any prefactors for the transormed argument u
 *
 * @param [in]  u          Transformed argument. Values in the interval [-1;1] are allowed
 * @param [out] integrands Logarithms and signs of the integrands
 * @param [in]  context    Context of the Drude-Integration
 */
static inline void integrands_drude_u(float80 u, integrands_drude_t *integrands,
                                      struct integ_context* context,
                                      unsigned int index,
                                      unsigned int iteration)
{
    if( fabs80(u - 1.0) < LDBL_EPSILON )
    {
        integrands->lnA_TE = -INFINITY;
        integrands->lnA_TM = -INFINITY;
        integrands->lnB_TE = -INFINITY;
        integrands->lnB_TM = -INFINITY;
        integrands->lnC_TE = -INFINITY;
        integrands->lnC_TM = -INFINITY;
        integrands->lnD_TE = -INFINITY;
        integrands->lnD_TM = -INFINITY;

        integrands->sign_A = +1;
        integrands->sign_B = +1;
        integrands->sign_C = +1;
        integrands->sign_D = +1;
    }
    else
    {
        const float80 x         = context->c0 * pow80_alpha((1+u) / (1-u));
        const float80 jacobi    = (pow80_alpha_m(1+u) / pow80_alpha_p(1-u) * (2 * ALPHA * context->c0));
        const float80 ln_jacobi = log80(fabs80(jacobi));

        integrands_drude(x, integrands, context, index, iteration);
    
        integrands->lnA_TM += ln_jacobi;
        integrands->lnA_TE += ln_jacobi;

        integrands->lnB_TM += ln_jacobi;
        integrands->lnB_TE += ln_jacobi;

        integrands->lnC_TM += ln_jacobi;
        integrands->lnC_TE += ln_jacobi;

        integrands->lnD_TM += ln_jacobi;
        integrands->lnD_TE += ln_jacobi;
    }
}


/*
 * Inititialize an object of the type struct drude_error.
 */
static inline void init_accuracy(struct drude_error* error,
                                 float80 ln_prefactor_A,
                                 float80 ln_prefactor_B,
                                 float80 ln_prefactor_CD)
{
    /*
     * If one prefactor is infinity (or -infinity) the iteration will stop immediately.
     * It doesn't matter. At the end the result will be infinity (or -infinity).
     */
    if(isinf(ln_prefactor_A))
        error->lnA = INFINITY;
    else
        error->lnA = log80(ACCURACY);

    if(isinf(ln_prefactor_B))
        error->lnB = INFINITY;
    else
        error->lnB = log80(ACCURACY);

    if(isinf(ln_prefactor_CD))
        error->lnC = error->lnD = INFINITY;
    else
    {
        error->lnC = log80(ACCURACY);
        error->lnD = log80(ACCURACY);
    }
}


/*
 * Set the Value of A,B,C,D in an integrands_drude_t object to infinity
 */
static inline void init_infty(integrands_drude_t* id)
{
    id->lnA_TE = id->lnA_TM = id->lnB_TE = id->lnB_TM = INFINITY;
    id->lnC_TE = id->lnC_TM = id->lnD_TE = id->lnD_TM = INFINITY;
    id->sign_A = id->sign_B = id->sign_C = id->sign_D = 1;
}


/*
 * Set the Value of A,B,C,D in an integrands_drude_t object to zero
 */
static inline void init_zero(integrands_drude_t* id)
{
    id->lnA_TE = id->lnA_TM = id->lnB_TE = id->lnB_TM = -INFINITY;
    id->lnC_TE = id->lnC_TM = id->lnD_TE = id->lnD_TM = -INFINITY;
    id->sign_A = id->sign_B = id->sign_C = id->sign_D = 1;
}


/*
 * Add Values of a integrands_drude_t to another one.
 * Something like "left += right" for integrands_drude_t.
 */
static inline void drude_plusequal(integrands_drude_t* left,
                                   integrands_drude_t* right)
{
    left->lnA_TE = logadd_s(left->lnA_TE, left->sign_A, right->lnA_TE,
                            right->sign_A, &left->sign_A);

    left->lnB_TE = logadd_s(left->lnB_TE, left->sign_B, right->lnB_TE,
                            right->sign_B, &left->sign_B);

    left->lnC_TE = logadd_s(left->lnC_TE, left->sign_C, right->lnC_TE,
                            right->sign_C, &left->sign_C);

    left->lnD_TE = logadd_s(left->lnD_TE, left->sign_D, right->lnD_TE,
                            right->sign_D, &left->sign_D);

    left->lnA_TM = logadd_s(left->lnA_TM, left->sign_A, right->lnA_TM,
                            right->sign_A, &left->sign_A);

    left->lnB_TM = logadd_s(left->lnB_TM, left->sign_B, right->lnB_TM,
                            right->sign_B, &left->sign_B);

    left->lnC_TM = logadd_s(left->lnC_TM, left->sign_C, right->lnC_TM,
                            right->sign_C, &left->sign_C);

    left->lnD_TM = logadd_s(left->lnD_TM, left->sign_D, right->lnD_TM,
                            right->sign_D, &left->sign_D);
}


/*
 * Add the Values of two integrands_drude_t.
 * Something like "result = left_arg + right_arg" for integrands_drude_t
 */
static inline void drude_plus(integrands_drude_t* result,
                              integrands_drude_t* left_arg,
                              integrands_drude_t* right_arg)
{
    result->lnA_TE = logadd_s(left_arg->lnA_TE, left_arg->sign_A, right_arg->lnA_TE,
                            right_arg->sign_A, &result->sign_A);

    result->lnB_TE = logadd_s(left_arg->lnB_TE, left_arg->sign_B, right_arg->lnB_TE,
                            right_arg->sign_B, &result->sign_B);

    result->lnC_TE = logadd_s(left_arg->lnC_TE, left_arg->sign_C, right_arg->lnC_TE,
                            right_arg->sign_C, &result->sign_C);

    result->lnD_TE = logadd_s(left_arg->lnD_TE, left_arg->sign_D, right_arg->lnD_TE,
                            right_arg->sign_D, &result->sign_D);


    result->lnA_TM = logadd_s(left_arg->lnA_TM, left_arg->sign_A, right_arg->lnA_TM,
                              right_arg->sign_A, &result->sign_A);

    result->lnB_TM = logadd_s(left_arg->lnB_TM, left_arg->sign_B, right_arg->lnB_TM,
                            right_arg->sign_B, &result->sign_B);

    result->lnC_TM = logadd_s(left_arg->lnC_TM, left_arg->sign_C, right_arg->lnC_TM,
                            right_arg->sign_C, &result->sign_C);

    result->lnD_TM = logadd_s(left_arg->lnD_TM, left_arg->sign_D, right_arg->lnD_TM,
                            right_arg->sign_D, &result->sign_D);
}


static inline void drude_minus(integrands_drude_t* result,
                               integrands_drude_t* left_arg,
                               integrands_drude_t* right_arg)
{
    /*
     * Subtract the Values of two integrands_drude_t.
     * Something like "result = left_arg - right_arg" for integrands_drude_t
     */
    right_arg->sign_A *= -1;
    right_arg->sign_B *= -1;
    right_arg->sign_C *= -1;
    right_arg->sign_D *= -1;
    drude_plus(result, left_arg, right_arg);
}


/*
 * Multiply the Values of an integrands_drude_t object by a (positive) factor.
 * "log_factor" is the logarithm of the factor.
 * id = factor * id
 */
static inline void drude_mult_factor(integrands_drude_t* id, float80 log_factor)
{
    id->lnA_TE += log_factor;
    id->lnA_TM += log_factor;

    id->lnB_TE += log_factor;
    id->lnB_TM += log_factor;

    id->lnC_TE += log_factor;
    id->lnC_TM += log_factor;

    id->lnD_TE += log_factor;
    id->lnD_TM += log_factor;
}


/** @brief Evaluate the integrals for the simple- and romberg-method
 *
 * This Function calculates the integrals A,B,C,D for the simple- and
 * romberg-method with all the necessary prefactors.
 *
 * @param [in]  self       Casimir object
 * @param [out] cint       Logarithms of values and signs of integrals
 * @param [in]  l1         Value of l1
 * @param [in]  l2         Value of l2
 * @param [in]  m          Value of m
 * @param [in]  n          Index of the Matsubara-Term
 * @param [in]  T          Scaled Temperature
 */
static inline void do_integrate(casimir_t* self, casimir_integrals_t* cint,
                                int l1, int l2, int m, int n, double T)
{
//    const float80 c0  = (l1 + l2 - 1) / (n*T);
    const float80 c0  = l1 + l2 - 1;
    const float80 tau = 2 * n * T;
    float80 c_max;
    integrands_drude_t total;

    const float80 ln_lambda = casimir_lnLambda(l1, l2, m, NULL); /* sign: -1 */
    const float80 log_m   = log80(m);
    const float80 log_tau = log80(tau);

    struct drude_error error;

    if(fabs80(c0 - 1.0) >= LDBL_EPSILON)
    {
        c_max = pow80( (10.0*c0 / (c0 - 1.0)), 1.0 / ALPHA);
        c_max = (c_max - 1) / (c_max + 1);
    }
    else
        c_max = 1.0;

    struct integ_context context = {
        .l1      = l1,
        .l2      = l2,
        .m       = m,
        .nT      = n*T,
        .c0      = c0,
        .c_max   = c_max,
        .casimir = self
    };


    /*
     * Integral C and D have the same prefactor.
     */
    const float80 prefactor_A     = ln_lambda + 2 * log_m  + log_tau - tau;
    const float80 prefactor_B     = ln_lambda - tau - 3 * log_tau;
    const float80 prefactor_CD    = ln_lambda + log_m - tau - log_tau;

    /*
     * First, we calculate the accuracy which is required
     * If one of the prefactors is -infinity (e.g. when m = 0), we don't
     * need to check if the integral converges. The result will be -infinity.
     */
    init_accuracy(&error, prefactor_A, prefactor_B, prefactor_CD);

    /*
     * First, the plm_cache needs to be initialized. When we change l1,l2,m or nT
     * we have to call the free function.
     */
#ifdef USE_PLM_CACHE
    plm_cache_init(&context, n);
#endif

    integrate_romberg(&context, &total, &error);

#ifdef USE_PLM_CACHE
    plm_cache_free(&context);
#endif


    // ---------------------------------------

    /* A */
    cint->lnA_TE = prefactor_A + total.lnA_TE;
    cint->lnA_TM = prefactor_A + total.lnA_TM;

    /* r_TE is negative, r_TM is positive and Lambda(l1,l2,m) is negative.
       => TM negative sign, TE positive sign */
    cint->signA_TM = -MPOW(l2+m) * total.sign_A;
    cint->signA_TE = +MPOW(l2+m) * total.sign_A;

    /* B */
    cint->lnB_TE = prefactor_B + total.lnB_TE;
    cint->lnB_TM = prefactor_B + total.lnB_TM;

    cint->signB_TM = -MPOW(l2+m+1) * total.sign_B;
    cint->signB_TE = +MPOW(l2+m+1) * total.sign_B;

    /* C */
    cint->lnC_TE = prefactor_CD + total.lnC_TE;
    cint->lnC_TM = prefactor_CD + total.lnC_TM;

    cint->signC_TM = -MPOW(l2+m) * total.sign_C;
    cint->signC_TE = +MPOW(l2+m) * total.sign_C;

    /* D */
    /* Same prefactor as C */
    cint->lnD_TE = prefactor_CD + total.lnD_TE;
    cint->lnD_TM = prefactor_CD + total.lnD_TM;

    cint->signD_TM = -MPOW(l2+m+1) * total.sign_D;
    cint->signD_TE = +MPOW(l2+m+1) * total.sign_D;
}


/*
 * ******************** Romberg Integration ********************
 */

/*
 * Use Romberg's method to integrate.
 * The Value of the Integral is I = lim_{k \rarrow \infty} I_{1,k}
 * In each Iteration, we will calculate a new row of the Terms I_{n,k}
 *
 * The approach is:
 * 1st Iteration    I11
 * 2nd Iteration    I21    I12
 * 3rd Iteration    I31    I22    I13
 * 4th Iteration    I41    I32    I23    I14
 * 5th Iteration    I51    I42    I33    I24    I15
 *     ...          ...    ...    ...    ...    ...
 */

/*
 * The romberg_chain we use is: 1, 1/2, 1/4, 1/8, ...
 * The Index starts with 1
 */
static inline float80 romberg_chain(unsigned int index) {
    return pow80(0.5,index-1);
}


/*
 * Calculate the Term I_{1,1}
 */
static void romberg_I11(integrands_drude_t* result,
                        struct integ_context* context)
{
    /*
     * The Term I_{1,1} is \frac{h_1}{2}(f(a) + f(b)).
     * In our case a = 0 and f(a)=0, so we can neglect that.
     */
    integrands_drude_u(context->c_max, result, context, 0, 1);
    drude_mult_factor(result, log80(0.5 * (context->c_max - (-1.0))));
}


/*
 * Calculate the Term I_{n,1}
 * Summarized Trapezoidal-rule with several intervals
 *
 * This function calculates:
 * I_{n,1} = \frac{h_n}{2} \left( f(a) + f(b) + 2 sum_{i=1}^{(h_1 / h_n) - 1} f(a + i h_n) \right)
 */
static inline void romberg_In1(integrands_drude_t* I_current,
                               integrands_drude_t* I_last,
                               struct integ_context* context,
                               unsigned int n)
{
    integrands_drude_t temp;
    float80 log_hn = log80(romberg_chain(n) * (context->c_max - (-1.0)));

    /*
     * Divide the Interval into steps of the size 1 / (2^n)
     * We have got 2^(n-2) points
     */
    const size_t num_points = pow(2, n-2);
    const float80 stepsize  = pow(0.5, n-1); 

    /*
     * We can reuse the points from the last iteration.
     * "I_last" is the sum evaluated for 2^(n-1) points.
     * We reuse theese values and add only the necessary points
     */
    *I_current = *I_last;
    
    drude_mult_factor(I_current, log80(0.5));

    for(unsigned int i = 1; i <= num_points; ++i)
    {
        /*
         * We already calculated the even points in the last iteration
         * We need only the odd ones (1*stepsize, 3*stepsize, 5*stepsize ...)
         */
        const float80 frac = ((float80)(2*i - 1)) * stepsize;
        const float80 u = frac * (context->c_max - (-1.0)) - 1.0;
        integrands_drude_u(u, &temp, context, i, n);
        drude_mult_factor(&temp, log_hn);
        drude_plusequal(&I_current[0], &temp);
    }
}


/*
 * Calculate I_{n,k} from I_{n+1,k-1} and I_{n,k-1}
 *
 * We use the formula:
 * I_{n,k} = \frac{2^{2(k-1)} I_{n+1,k-1} - I_{n,k-1}}{2^{2(k-1)} - 1}
 */
static inline void romberg_subcycle(integrands_drude_t* a,
                                    integrands_drude_t* b,
                                    integrands_drude_t* c,
                                    unsigned int k)
{
    /*
     * a: I_{n,k}
     * b: I_{n+1,k-1}
     * c: I_{n,k-1}
     */
    const float80 prefactor = pow80(2, 2 * (k-1));
    
    *a = *b;
    drude_mult_factor(a, log80(prefactor));
    drude_minus(a, a, c);
    drude_mult_factor(a, log80(1 / (prefactor - 1.0)));
}


/*
 * Calculate the Accuracy of the Integration.
 * "cur_result" needs to be I_{1,n+1}.
 * "last_results" needs to be I_{1,n}.
 * Returns 1 if the result is more accurate than "error". Else return 0.
 */
static inline int romberg_is_accurate(integrands_drude_t* cur_result,
                                      integrands_drude_t* last_result,
                                      struct drude_error* error)
{
    float80 tmp;

    tmp = logdiff_rel(cur_result->lnA_TE, last_result->lnA_TE);
    if(tmp > error->lnA || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnA_TM, last_result->lnA_TM);
    if(tmp > error->lnA || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnB_TE, last_result->lnB_TE);
    if(tmp > error->lnB || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnB_TM, last_result->lnB_TM);
    if(tmp > error->lnB || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnC_TE, last_result->lnC_TE);
    if(tmp > error->lnC || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnC_TM, last_result->lnC_TM);
    if(tmp > error->lnC || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnD_TE, last_result->lnD_TE);
    if(tmp > error->lnD || isnan(tmp))
        return 0;

    tmp = logdiff_rel(cur_result->lnD_TM, last_result->lnD_TM);
    if(tmp > error->lnD || isnan(tmp))
        return 0;

    return 1;
}


/*
 * @brief Romberg integration function. Evalues the integrals without any prefactors.
 *
 * @param [in]  context   Drude Integration Context
 * @param [out] result    Value of the 8 integrals (A_TE, A_TM, B_TE, B_TM, ...)
 * @param [in]  error     Claimed accuracy
 */
static inline void integrate_romberg(struct integ_context* context,
                                        integrands_drude_t* result,
                                        struct drude_error* error)
{
    const size_t max_order = 20;
    const size_t min_order = 5;
    integrands_drude_t I1[max_order], I2[max_order];
    /*
     * I_current stores the terms In,k of this cycle, while I_last stores the
     * terms of the last cycle.
     * So e.g.: cur_cycle: I4,1  I3,2  I2,3  I1,4
     *         last_cycle: I3,1  I2,2  I1,3
     */
    integrands_drude_t* I_current, *I_last;
    unsigned int n, k;

    n = 1;
    k = 1;

    I_current = I1;
    I_last    = I2;

    romberg_I11(&I_last[0], context);
    for(n = 2; n <= max_order; ++n)
    {
        // Calculate I(n, 1)
        romberg_In1(I_current, I_last, context, n);
        // Calculate the Subcycle
        for(k = 2; k <= n; ++k)
        {
            romberg_subcycle(&I_current[k - 1], &I_current[k - 2], &I_last[k - 2], k);
        }

        /*
         * Swap I_current and I_last.
         * I_current is now the array of the last main cycle.
         * I_last is the array for the new cycle.
         */
        integrands_drude_t* id_temp;
        id_temp   = I_current;
        I_current = I_last;
        I_last    = id_temp;

        /*
         * Check if the result is accurate enough and exit in this case.
         * Do at least min_order Iterations.
         * The error can be calculated with I_{1,n+1} and I_{1,n}.
         * Keep in mind: we already swapped I_current and I_last.
         */
        if(n >= min_order &&
           romberg_is_accurate(&I_last[n-1], &I_current[n-2], error))
        {
            n++;
            
            /*
             * Integral converged. Return the result of the last iteration.
             */
            *result = I_last[n - 2];
            return;
        }
    }
    /*
     * We did "max_order" iterations and the integral is still not accurate
     * enough.
     */
    WARN(1, "Drude-Integration: Failed to achieve accuracy of %g.", ACCURACY);

    /*
     * Leave this function and the return the result of the last iteration.
     */
    *result = I_last[n - 2];
}

#endif



/** @brief Calculate integrals A,B,C,D including prefactor Lambda vor Drude metals
 *
 * This function calculates
 *    Lambda(l1,l2,m)*A_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*B_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*C_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*D_(l1,l2)^(m)
 * for Drude metals.
 *
 * @param [in]  self Casimir object
 * @param [out] cint logarithms of values and signs of integrals
 * @param [in]  l1   \f$\ell_1\f$
 * @param [in]  l2   \f$\ell_2\f$
 * @param [in]  m    \f$m\f$
 * @param [in]  nT   \f$nT\f$
 */
void casimir_integrate_drude(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, int n, double T)
{
#if defined INTEGRATION_GAUSS_LAGUERRE
    integrate_gauss_laguerre(self, cint, l1, l2, m, n, T);
#else
    do_integrate(self, cint, l1, l2, m, n, T);
#endif
}


void casimir_integrate_drude_init(casimir_t* self)
{
#if (defined INTEGRATION_ROMBERG && defined USE_PLM_CACHE)
    plm_create_cache(self);
#endif
}


void casimir_integrate_drude_free()
{
#if (defined INTEGRATION_ROMBERG && defined USE_PLM_CACHE)
    plm_destroy_cache();
#endif    
}


/*
 * ******************** Gauss-Laguerre Integration ********************
 */
#ifdef INTEGRATION_GAUSS_LAGUERRE

static void integrands_drude_laguerre(float80 x, integrands_drude_t *integrands,
                                      casimir_t *self, double nT, int l1, int l2, int m)
{
    plm_combination_t comb;
    const float80 tau        = 2*nT;
    const float80 k          = sqrt80(pow_2(x)/4 + nT*x);
    const float80 log_factor = log80(pow_2(x)+2*tau*x);
    float80 r_TE, r_TM;

    casimir_rp(self, nT, k, &r_TE, &r_TM);
    const float80 lnr_TE = log80(-r_TE);
    const float80 lnr_TM = log80(r_TM);

    plm_PlmPlm(l1, l2, m, 1+x/tau, &comb);

    const float80 A = comb.lnPl1mPl2m - log_factor;
    integrands->lnA_TE = lnr_TE + A;
    integrands->lnA_TM = lnr_TM + A;
    integrands->sign_A = comb.sign_Pl1mPl2m;

    const float80 B = comb.lndPl1mdPl2m + log_factor;
    integrands->lnB_TE = lnr_TE + B;
    integrands->lnB_TM = lnr_TM + B;
    integrands->sign_B = comb.sign_dPl1mdPl2m;

    const float80 C = comb.lnPl1mdPl2m;
    integrands->lnC_TE = lnr_TE + C;
    integrands->lnC_TM = lnr_TM + C;
    integrands->sign_C = comb.sign_Pl1mdPl2m;

    const float80 D = comb.lndPl1mPl2m;
    integrands->lnD_TE = lnr_TE + D;
    integrands->lnD_TM = lnr_TM + D;
    integrands->sign_D = comb.sign_dPl1mPl2m;
}


/** @brief Calculate integrals A,B,C,D including prefactor Lambda vor Drude metals
 *
 * This function calculates
 *    Lambda(l1,l2,m)*A_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*B_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*C_(l1,l2)^(m),
 *    Lambda(l1,l2,m)*D_(l1,l2)^(m)
 * for Drude metals.
 *
 * @param [in]  self Casimir object
 * @param [out] cint logarithms of values and signs of integrals
 * @param [in]  l1   \f$\ell_1\f$
 * @param [in]  l2   \f$\ell_2\f$
 * @param [in]  m    \f$m\f$
 * @param [in]  nT   \f$nT\f$
 */
void integrate_gauss_laguerre(casimir_t *self, casimir_integrals_t *cint, int l1, int l2, int m, int n, double T)
{
    integrands_drude_t integrand;
    const float80 tau = 2*n*T;
    const float80 ln_tau = log80(tau);
    const float80 ln_Lambda = casimir_lnLambda(l1, l2, m, NULL); /* sign: -1 */
    float80 prefactor;
    float80 *xk, *ln_wk;
    /* XXX use order 50 at the moment XXX */
    const int N = gausslaguerre_nodes_weights(50, &xk, &ln_wk);

    /* allocate space for signs_A, signs_B, signs_C, signs_D */
    sign_t *signs_ABCD = xmalloc(4*N*sizeof(sign_t));
    sign_t *signs_A = signs_ABCD;
    sign_t *signs_B = signs_A+1*N;
    sign_t *signs_C = signs_A+2*N;
    sign_t *signs_D = signs_A+3*N;

    /* allocate space for lnA_TE, lnA_TM, lnB_TE, lnB_TM, lnC_TE, lnC_TM,
     * lnD_TE, lnD_TM */
    float80 *ln_ABCD = xmalloc(4*2*N*sizeof(integrands_drude_t));
    float80 *lnA_TE  = ln_ABCD;
    float80 *lnA_TM  = ln_ABCD + 1*N;
    float80 *lnB_TE  = ln_ABCD + 2*N;
    float80 *lnB_TM  = ln_ABCD + 3*N;
    float80 *lnC_TE  = ln_ABCD + 4*N;
    float80 *lnC_TM  = ln_ABCD + 5*N;
    float80 *lnD_TE  = ln_ABCD + 6*N;
    float80 *lnD_TM  = ln_ABCD + 7*N;

    for(int i = 0; i < N; i++)
    {
        integrands_drude_laguerre(xk[i], &integrand, self, n*T, l1, l2, m);

        lnA_TE[i]  = ln_wk[i] + integrand.lnA_TE;
        lnA_TM[i]  = ln_wk[i] + integrand.lnA_TM;
        signs_A[i] = integrand.sign_A;

        lnB_TE[i]  = ln_wk[i] + integrand.lnB_TE;
        lnB_TM[i]  = ln_wk[i] + integrand.lnB_TM;
        signs_B[i] = integrand.sign_B;

        lnC_TE[i]  = ln_wk[i] + integrand.lnC_TE;
        lnC_TM[i]  = ln_wk[i] + integrand.lnC_TM;
        signs_C[i] = integrand.sign_C;

        lnD_TE[i]  = ln_wk[i] + integrand.lnD_TE;
        lnD_TM[i]  = ln_wk[i] + integrand.lnD_TM;
        signs_D[i] = integrand.sign_D;
    }


    /* B */
    prefactor = ln_Lambda -tau-3*ln_tau; /* exp(-tau)/tau³ */
    cint->lnB_TE = prefactor + logadd_ms(lnB_TE, signs_B, N, &cint->signB_TE);
    cint->lnB_TM = prefactor + logadd_ms(lnB_TM, signs_B, N, &cint->signB_TM);

    cint->signB_TM = -MPOW(l2+m+1) * cint->signB_TM;
    cint->signB_TE = +MPOW(l2+m+1) * cint->signB_TE;


    if(m > 0)
    {
        const float80 log_m = log(m);

        /* A */
        prefactor = ln_Lambda + 2*log_m+ln_tau-tau; /* m²*tau*exp(-tau) */
        cint->lnA_TE = prefactor + logadd_ms(lnA_TE, signs_A, N, &cint->signA_TE);
        cint->lnA_TM = prefactor + logadd_ms(lnA_TM, signs_A, N, &cint->signA_TM);

        /* r_TE is negative, r_TM is positive and Lambda(l1,l2,m) is negative.
           => TM negative sign, TE positive sign */
        cint->signA_TM = -MPOW(l2+m) * cint->signA_TM;
        cint->signA_TE = +MPOW(l2+m) * cint->signA_TE;


        /* C */
        prefactor = ln_Lambda + log_m-tau-ln_tau; /* m*exp(-tau)/tau */
        cint->lnC_TE = prefactor + logadd_ms(lnC_TE, signs_C, N, &cint->signC_TE);
        cint->lnC_TM = prefactor + logadd_ms(lnC_TM, signs_C, N, &cint->signC_TM);

        cint->signC_TM = -MPOW(l2+m) * cint->signC_TM;
        cint->signC_TE = +MPOW(l2+m) * cint->signC_TE;


        /* D */
        /* prefactor is identical to C */
        cint->lnD_TE = prefactor + logadd_ms(lnD_TE, signs_D, N, &cint->signD_TE);
        cint->lnD_TM = prefactor + logadd_ms(lnD_TM, signs_D, N, &cint->signD_TM);

        cint->signD_TM = -MPOW(l2+m+1) * cint->signD_TM;
        cint->signD_TE = +MPOW(l2+m+1) * cint->signD_TE;
    }
    else
    {
        cint->lnA_TM = cint->lnA_TE = -INFINITY;
        cint->signA_TM = cint->signA_TE = +1;

        cint->lnC_TM = cint->lnC_TE = -INFINITY;
        cint->signC_TM = cint->signC_TE = +1;

        cint->lnD_TM = cint->lnD_TE = -INFINITY;
        cint->signD_TM = cint->signD_TE = +1;
    }

    xfree(ln_ABCD);
    xfree(signs_ABCD);
}
#endif
