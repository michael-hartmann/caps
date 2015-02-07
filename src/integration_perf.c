#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_perf.h"

/* p must have length len_p1+len_p2-1 */
static void polymult(edouble p1[], int len_p1, edouble p2[], int len_p2, edouble p[])
{
    int i,j;

    for(i = 0; i < len_p1+len_p2-1; i++)
        p[i] = 0;

    for(i = 0; i < len_p1; i++)
        for(j = 0; j < len_p2; j++)
            p[i+j] += p1[i]*p2[j];
}


/* p must have length m+n+1 */
static void poly1(int m, int n, edouble p[])
{
    /* (z+2)^m * (z+1)^n */
    edouble p1[m+1];
    edouble p2[n+1];
    int k;

    for(k = 0; k < m+1; k++)
        p1[k] = 0;
    for(k = 0; k < n+1; k++)
        p2[k] = 0;

    for(k = 0; k <= m; k++)
        p1[k] = expq(lgammaq(m+1)-lgammaq(k+1)-lgammaq(m+1-k)+(m-k)*LOG2);

    for(k = 0; k <= n; k++)
        p2[k] = expq(lgammaq(n+1)-lgammaq(k+1)-lgammaq(n+1-k));

    polymult(p1,m+1,p2,n+1,p);
}


/* p must have size nu+1-2m */
static void poly2(int nu, int m2, edouble p[])
{
    /* m2 = 2*m
     * d^2m/dz^2m P_(nu)(1+z)
     */
    //const int len = nu-m2+1;
    int k;

    for(k = m2; k <= nu; k++)
        p[k-m2] = expq(lgammaq(k+nu+1)-lgammaq(k+1)-lgammaq(k-m2+1)-lgammaq(-k+nu+1)-k*LOG2);
}


static edouble polyintegrate(edouble p[], int len_p, int offset, edouble tau)
{
    int k;
    edouble value = 0;

    // XXX prevent overflows
    for(k = offset; k < len_p+offset; k++)
        //value += gamma(k+1)/tau**(k+1)*p[k-offset]
        value += expq(lgammaq(k+1)-(k+1)*logq(tau))*p[k-offset];

    return value;
}


static edouble I(int beta, int nu, int m2, edouble tau)
{
    // exp(-z*tau) * (z^2+2z)^-1 * (1+z)^beta * Plm(nu, 2m, 1+z)
    int m = m2/2;
    edouble p1[m+beta], p2[nu+1-m2], p[-m+beta+nu];

    poly1(m-1, beta, p1);
    poly2(nu,m2,p2);
    polymult(p1, m+beta, p2, nu+1-m2, p);

    return polyintegrate(p, -m+beta+nu, m-1, tau);
}

static void polydplm(edouble pl1[], edouble pl2[], int l1, int l2, int m)
{
    int k;

    if(m == 0)
    {
        for(k = 1; k <= l1; k++)
            pl1[k-1] = expq(lngamma(1+k+l1)-lngamma(1+k)-lngamma(1+l1-k)-lngamma(k)-(k-1)*LOG2)/2;

        for(k = 1; k <= l2; k++)
            pl2[k-1] = expq(lngamma(1+k+l2)-lngamma(1+k)-lngamma(1+l2-k)-lngamma(k)-(k-1)*LOG2)/2;
    }
    else
    {
        for(k = 0; k <= l1-m+1; k++)
            pl1[k] = 0;
        for(k = 0; k <= l2-m+1; k++)
            pl2[k] = 0;

        pl1[l1+1-m] += (l1-m+1)*expq(lngamma(2*l1+3)-lngamma(l1+2)-lngamma(l1-m+2)-(l1+1-m)*LOG2);
        for(k = 0; k <= l1-m; k++)
        {
            edouble common = expq(lngamma(1+k+l1+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l1-k-m)-k*LOG2);
            pl1[k] += common*(pow_2(m)+m*(k-l1-1)-2.0*k*l1-2*k)/(m-l1+k-1);
            pl1[k+1] -= common*(l1+1);
        }

        pl2[l2+1-m] += (l2-m+1)*expq(lngamma(2*l2+3)-lngamma(l2+2)-lngamma(l2-m+2)-(l2+1-m)*LOG2);
        for(k = 0; k <= l2-m; k++)
        {
            edouble common = expq(lngamma(1+k+l2+m)-lngamma(1+k+m)-lngamma(1+k)-lngamma(1+l2-k-m)-k*LOG2);
            pl2[k] += common*(pow_2(m)+m*(k-l2-1)-2.0*k*l2-2*k)/(m-l2+k-1);
            pl2[k+1] -= common*(l2+1.);
        }
    }
}


/* TODO: m = 0, check signs! */
void casimir_integrate_perf(casimir_integrals_t *cint, int l1, int l2, int m, double nT)
{
    int nu,q;
    edouble tau = 2*nT;
    edouble log_m;
    const int qmax_l1l2   = GAUNT_QMAX(l1,  l2,  m,m);
    edouble a_l1l2[qmax_l1l2+1];
    const int qmax_l1pl2  = GAUNT_QMAX(l1+1,l2,  m,m);
    edouble a_l1pl2[qmax_l1pl2+1];
    const int qmax_l1l2p  = GAUNT_QMAX(l1,  l2+1,m,m);
    edouble a_l1l2p[qmax_l1l2p+1];
    const int qmax_l1pl2p = GAUNT_QMAX(l1+1,l2+1,m,m);
    edouble a_l1pl2p[qmax_l1pl2p+1];
    edouble log_A,log_B,log_C,log_D;
    edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);
    int sign_A, sign_B, sign_C, sign_D;

    if(m == 0)
    {
        edouble value;
        edouble pm[] = {0,2,1}; /* z²+2z */
        edouble pdpl1m[l1-m+2];
        edouble pdpl2m[l2-m+2];
        edouble pdpl1mpdpl2m[l1+l2-1];
        edouble pmpdpl1mpdpl2m[l1+l2+1];

        polydplm(pdpl1m,pdpl2m,l1,l2,m);

        polymult(pdpl1m, l1, pdpl2m, l2, pdpl1mpdpl2m);
        polymult(pm, 3, pdpl1mpdpl2m, l1+l2-1, pmpdpl1mpdpl2m);

        value = polyintegrate(pmpdpl1mpdpl2m, l1+l2+1, 0, tau);

        cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+logq(fabsq(value));
        cint->signB_TM = -MPOW(l2+1)*copysignq(1,value);
        cint->signB_TE = -cint->signB_TM;

        cint->lnA_TM = cint->lnA_TE = -INFINITY;
        cint->lnC_TM = cint->lnC_TE = -INFINITY;
        cint->lnD_TM = cint->lnD_TE = -INFINITY;
        cint->signA_TM = cint->signA_TE = +1;
        cint->signC_TM = cint->signC_TE = +1;
        cint->signD_TM = cint->signD_TE = +1;

        return;
    }

    /* calculate Gaunt coefficients */
    gaunt(l1,   l2,   m, m, a_l1l2);
    gaunt(l1+1, l2,   m, m, a_l1pl2);
    gaunt(l1,   l2+1, m, m, a_l1l2p);
    gaunt(l1+1, l2+1, m, m, a_l1pl2p);

    /* A */
    {
        edouble A = 0;
        for(q = 0; q <= qmax_l1l2; q++)
        {
            nu = l1+l2-2*q;
            A += a_l1l2[q]*I(0,nu,2*m,tau);
        }
        A *= gaunt_a0(l1,l2,m,m);

        log_A  = logq(fabsq(A));
        sign_A = copysignq(1, A);
    }

    /* B */
    {
        edouble B, B1, B2, B3, B4;

        B1 = 0;
        for(q = 0; q <= qmax_l1pl2p; q++)
        {
            nu = l1+1+l2+1-2*q;
            B1 += a_l1pl2p[q]*I(0,nu,2*m,tau);
        }
        B1 *= gaunt_a0(l1+1,l2+1,m,m);
        B1 *= (l1-m+1)*(l2-m+1);

        B2 = 0;
        for(q = 0; q <= qmax_l1pl2; q++)
        {
            nu = l1+1+l2-2*q;
            B2 += a_l1pl2[q]*I(1,nu,2*m,tau);
        }
        B2 *= -gaunt_a0(l1+1,l2,m,m);
        B2 *= (l1-m+1)*(l2+1);

        B3 = 0;
        for(q = 0; q <= qmax_l1l2p; q++)
        {
            nu = l1+l2+1-2*q;
            B3 += a_l1l2p[q]*I(1,nu,2*m,tau);
        }
        B3 *= -gaunt_a0(l1,l2+1,m,m);
        B3 *= (l1+1)*(l2-m+1);

        B4 = 0;
        for(q = 0; q <= qmax_l1l2; q++)
        {
            nu = l1+l2-2*q;
            B4 += a_l1l2[q]*I(2,nu,2*m,tau);
        }
        B4 *= gaunt_a0(l1,l2,m,m);
        B4 *= (l1+1)*(l2+1);

        B = B1+B2+B3+B4;

        log_B  = logq(fabsq(B));
        sign_B = copysignq(1,B);
    }

    /* C */
    {
        edouble C,C1,C2;

        C1 = 0;
        for(q = 0; q <= qmax_l1l2p; q++)
        {
            nu = l1+l2+1-2*q;
            C1 += a_l1l2p[q]*I(0,nu,2*m,tau);
        }
        C1 *= gaunt_a0(l1,l2+1,m,m)*(l2-m+1);

        C2 = 0;
        for(q = 0; q <= qmax_l1l2; q++)
        {
            nu = l1+l2-2*q;
            C2 += a_l1l2[q]*I(1,nu,2*m,tau);
        }
        C2 *= -gaunt_a0(l1,l2,m,m)*(l2+1);

        C = C1+C2;
        log_C  = logq(fabsq(C));
        sign_C = copysignq(1,C);
    }

    /* D */
    {
        edouble D,D1,D2;

        D1 = 0;
        for(q = 0; q <= qmax_l1pl2; q++)
        {
            nu = l1+1+l2-2*q;
            D1 += a_l1pl2[q]*I(0,nu,2*m,tau);
        }
        D1 *= gaunt_a0(l1+1,l2,m,m)*(l1-m+1);

        D2 = 0;
        for(q = 0; q <= qmax_l1l2; q++)
        {
            nu = l1+l2-2*q;
            D2 += a_l1l2[q]*I(1,nu,2*m,tau);
        }
        D2 *= -gaunt_a0(l1,l2,m,m)*(l1+1);

        D = D1+D2;
        log_D  = logq(fabsq(D));
        sign_D = copysignq(1,D);
    }


    log_m   = logq(m);

    cint->lnA_TM   = cint->lnA_TE = 2*log_m+lnLambda-tau+log_A;
    cint->signA_TM = -MPOW(l2)*sign_A; /* - because Lambda(l1,l2,m) is negative */
    cint->signA_TE = -cint->signA_TM;

    cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+log_B;
    cint->signB_TM = -MPOW(l2+1)*sign_B;
    cint->signB_TE = -cint->signB_TM;

    cint->lnC_TM   = cint->lnC_TE = log_m+lnLambda-tau+log_C;
    cint->signC_TM = -MPOW(l2+1)*sign_C;
    cint->signC_TE = -cint->signC_TM;

    cint->lnD_TM   = cint->lnD_TE = log_m+lnLambda-tau+log_D;
    cint->signD_TM = -MPOW(l2)*sign_D; // XXX ???
    cint->signD_TE = -cint->signD_TM;
}
