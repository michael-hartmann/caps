#include <math.h>

#include "sfunc.h"
#include "gaunt.h"

/**
 * @brief Determine qmax
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$q_\mathrm{max}\f$
 */
int gaunt_qmax(const int n, const int nu, const int m)
{
    int xi = (n+nu-2*m)/2;
    return MIN(MIN(n,nu), xi);
}

/**
 * @brief Calculate \f$\log a_0\f$
 *
 * Cf. eq. (20).
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$\log a_0\f$
 */
double gaunt_log_a0(int n, int nu, int m)
{
    return lfac(2*n)-lfac(n)+lfac(2*nu)-lfac(nu)+lfac(n+nu)-lfac(2*n+2*nu)+lfac(n+nu-2*m)-lfac(n-m)-lfac(nu-m);
}

/**
 * @brief Calculate \f$a_0\f$
 *
 * Cf. eq. (20).
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @return a0 \f$a_0\f$
 */
double gaunt_a0(int n, int nu, int m)
{
    return exp(gaunt_log_a0(n,nu,m));
}

/* eq. (3) */
#define alpha(p, n, nu) (((pow_2(p)-pow_2(n+nu+1))*(pow_2(p)-pow_2(n-nu)))/(4*pow_2(p)-1))


/**
 * @brief Calculate Gaunt coefficients
 *
 * Determine Gaunt coefficients \f$a(m, n, mu, nu, p)\f$ for \f$m\f$, \f$n\f$,
 * \f$\mu\f$ and \f$\nu\f$ fixed.  These coefficients can be used to express
 * the product of two associated Legendre polynomials:
 *
 * \f[
 * P_n^m(x) P_{\nu}^{\mu}(x) = a_0 \sum_{q=0}^{q_\mathrm{max}} \tilde a_q P_{n+\nu-2q}^{m+mu}(x)
 * \f]
 *
 * \f$q_\mathrm{max}\f$ is the upper bound of summation, \f$a_0\f$ is the
 * prefactor and \f$\tilde a_q\f$ are normalized Gaunt coefficients.
 *
 * See [1] for more information, especially chapter 3. There is a brief
 * outline how to calculate Gaunt coefficients at the end of the chapter.
 *
 * Ref.: [1] Y.-L. Xu, J. Comp. Appl. Math. 85, 53 (1997)
 *
 * @param [in]  n  \f$n\f$
 * @param [in]  nu \f$\nu\f$
 * @param [in]  m  \f$m=\mu\f$
 * @param [out] a_tilde \f$\tilde a_q\f$ list of normalized Gaunt coefficients
 */
void gaunt(const int n_, const int nu_, const int m_, double a_tilde[])
{
    const double n  = n_;
    const double nu = nu_;
    const double m  = m_;
    const double n4 = n+nu-2*m;

    /* eq. (24) */
    const int qmax = gaunt_qmax(n,nu,m);

    /* eq. (28) */
    const double Ap = -2*m*(n-nu)*(n+nu+1);

    if(qmax < 0)
        return;

    a_tilde[0] = 1;
    if(qmax == 0)
        return;

    /* eq. (29) */
    a_tilde[1] = (n+nu-1.5)*(1-(2*n+2*nu-1)/(n4*(n4-1))*((m-n)*(m-n+1)/(2*n-1)+(m-nu)*(m-nu+1)/(2*nu-1)));
    if(qmax == 1)
        return;

    /* eq. (35) */
    a_tilde[2] = (2*n+2*nu-1)*(2*n+2*nu-7)/4*( (2*n+2*nu-3)/(n4*(n4-1)) * ( (2*n+2*nu-5)/(2*(n4-2)*(n4-3)) \
                * ( (m-n)*(m-n+1)*(m-n+2)*(m-n+3)/(2*n-1)/(2*n-3) \
                + 2*(m-n)*(m-n+1)*(m-nu)*(m-nu+1)/((2*n-1)*(2*nu-1)) \
                + (m-nu)*(m-nu+1)*(m-nu+2)*(m-nu+3)/(2*nu-1)/(2*nu-3) ) - (m-n)*(m-n+1)/(2*n-1) \
                - (m-nu)*(m-nu+1)/(2*nu-1) ) +0.5);

    for(int q = 3; q <= qmax; q++)
    {
        const int p = n+nu-2*q;
        const int p1 = p-2*m;
        const int p2 = p+2*m;

        if(Ap != 0)
        {
            /* eqs. (26), (27) */
            double c0 = (p+2)*(p+3)*(p1+1)*(p1+2)*Ap*alpha(p+1,n,nu);
            double c1 = Ap*(Ap*Ap \
               + (p+1)*(p+3)*(p1+2)*(p2+2)*alpha(p+2,n,nu) \
               + (p+2)*(p+4)*(p1+3)*(p2+3)*alpha(p+3,n,nu));
            double c2 = -(p+2)*(p+3)*(p2+3)*(p2+4)*Ap*alpha(p+4,n,nu);

            a_tilde[q] = (c1*a_tilde[q-1] + c2*a_tilde[q-2])/c0;
        }
        else
            /* eq. (30) */
            a_tilde[q] = (p+1)*(p2+2)*alpha(p+2,n,nu)*a_tilde[q-1] / ((p+2)*(p1+1)*alpha(p+1,n,nu));
    }
}
