#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_perf.h"

#include "utils.h"


/* p must have length m+n+1 */
static void poly1(int m, edouble p[])
{
    /* (z+2)^m */
    int k;

    for(k = 0; k <= m; k++)
        p[k] = expq(lgammaq(m+1)-lgammaq(k+1)-lgammaq(m+1-k)+(m-k)*LOG2);
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


/* evaluete integral I_nu^2m(tau) = exp(-z*tau)/ (z^2+2z) * Plm(nu, 2m, 1+z) */
static inline edouble I(integration_perf_t *self, int nu, int m2)
{
    const int m = m2/2;
    const int index = m*self->nu_max + nu;
    edouble v;

    v = self->cache_I[index];
    if(v == 0)
    {
        edouble tau = self->tau;

        edouble p1[m], p2[nu+1-m2], p[-m+nu];

        poly1(m-1, p1);
        poly2(nu,m2,p2);
        polymult(p1, m, p2, nu+1-m2, p);

        v = MPOW(m)*polyintegrate(p, -m+nu, m-1, tau);
        self->cache_I[index] = v;
    }

    return v;
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

void casimir_integrate_perf_init(integration_perf_t *self, double nT, int lmax)
{
    int i,size;
    self->tau = 2*nT;

    /* init cache */
    self->nu_max = 2*lmax+4;
    self->m2_max = 1*lmax+2;
    size = self->nu_max*self->m2_max;
    self->cache_I = xmalloc(size*sizeof(edouble));

    for(i = 0; i < size; i++)
        self->cache_I[i] = 0;

    self->cache_gaunt = xmalloc(sizeof(gaunt_cache_t));
    cache_gaunt_init(self->cache_gaunt, lmax);
}

void casimir_integrate_perf_free(integration_perf_t *self)
{
    cache_gaunt_free(self->cache_gaunt);
    xfree(self->cache_gaunt);

    xfree(self->cache_I);
    self->cache_I = NULL;
}

/* TODO: m = 0, check signs! */
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint, gaunt_cache_t *cache)
{
    edouble tau = self->tau;
    edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);

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
    }
    else
    {
        int nu,q, sign_A, sign_B, sign_C, sign_D;
        edouble log_m = logq(m);
        edouble log_A,log_B,log_C,log_D;

        const int qmax_l1l2   = gaunt_qmax(l1,  l2,  m);

        const int qmax_l1pl2  = gaunt_qmax(l1+1,l2,  m);
        const int qmax_l1ml2  = gaunt_qmax(l1-1,l2,  m);
        const int qmax_l1l2p  = gaunt_qmax(l1,  l2+1,m);
        const int qmax_l1l2m  = gaunt_qmax(l1,  l2-1,m);

        const int qmax_l1pl2p = gaunt_qmax(l1+1,l2+1,m);
        const int qmax_l1pl2m = gaunt_qmax(l1+1,l2-1,m);
        const int qmax_l1ml2p = gaunt_qmax(l1-1,l2+1,m);
        const int qmax_l1ml2m = gaunt_qmax(l1-1,l2-1,m);

        edouble *a_l1l2, *a_l1pl2, *a_l1ml2, *a_l1l2p, *a_l1l2m, *a_l1pl2p;
        edouble *a_l1pl2m, *a_l1ml2p, *a_l1ml2m;

        if(cache == NULL)
        {
            /* reserve space for gaunt coefficients on stack */
            a_l1l2 = alloca(MAX(1,qmax_l1l2+1)*sizeof(edouble));

            a_l1pl2 = alloca(MAX(1,qmax_l1pl2+1)*sizeof(edouble));
            a_l1ml2 = alloca(MAX(1,qmax_l1ml2+1)*sizeof(edouble));
            a_l1l2p = alloca(MAX(1,qmax_l1l2p+1)*sizeof(edouble));
            a_l1l2m = alloca(MAX(1,qmax_l1l2m+1)*sizeof(edouble));

            a_l1pl2p = alloca(MAX(1,qmax_l1pl2p+1)*sizeof(edouble));
            a_l1pl2m = alloca(MAX(1,qmax_l1pl2m+1)*sizeof(edouble));
            a_l1ml2p = alloca(MAX(1,qmax_l1ml2p+1)*sizeof(edouble));
            a_l1ml2m = alloca(MAX(1,qmax_l1ml2m+1)*sizeof(edouble));

            /* calculate Gaunt coefficients */
            gaunt(l1, l2,  m, a_l1l2);

            gaunt(l1-1, l2,  m, a_l1ml2);
            gaunt(l1,   l2-1,m, a_l1l2m);
            gaunt(l1-1, l2-1,m, a_l1ml2m);

            if((l1+1) >= m)
            {
                gaunt(l1+1, l2,  m, a_l1pl2);
                gaunt(l1+1, l2-1,m, a_l1pl2m);
            }

            if((l2+1) >= m)
            {
                gaunt(l1,   l2+1,m, a_l1l2p);
                gaunt(l1-1, l2+1,m, a_l1ml2p);
            }

            if((l1+1) >= m && (l2+1) >= m)
                gaunt(l1+1, l2+1, m, a_l1pl2p);
        }
        else
        {
            a_l1l2 = cache_gaunt_get(cache, l1,l2,m);

            a_l1pl2 = cache_gaunt_get(cache, l1+1,l2,m);
            a_l1ml2 = cache_gaunt_get(cache, l1-1,l2,m);
            a_l1l2p = cache_gaunt_get(cache, l1,l2+1,m);
            a_l1l2m = cache_gaunt_get(cache, l1,l2-1,m);

            a_l1pl2p = cache_gaunt_get(cache, l1+1,l2+1,m);
            a_l1pl2m = cache_gaunt_get(cache, l1+1,l2-1,m);
            a_l1ml2p = cache_gaunt_get(cache, l1-1,l2+1,m);
            a_l1ml2m = cache_gaunt_get(cache, l1-1,l2-1,m);
        }


        /* A */
        {
            edouble A = 0;
            for(q = 0; q <= qmax_l1l2; q++)
            {
                nu = l1+l2-2*q;
                A += a_l1l2[q]*I(self,nu,2*m);
            }

            log_A  = gaunt_log_a0(l1,l2,m)+logq(fabsq(A));
            sign_A = copysignq(1, A);
        }

        /* B */
        {
            edouble B, B1, B2, B3, B4;

            B1 = 0;
            if((l1-1) >= m && (l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2m; q++)
                {
                    nu = l1-1+l2-1-2*q;
                    B1 += a_l1ml2m[q]*I(self,nu,2*m);
                }
                B1 *= gaunt_a0(l1-1,l2-1,m);
                B1 *= (l1+1.)*(l1+m)*(l2+1.)*(l2+m);
            }

            B2 = 0;
            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2p; q++)
                {
                    nu = l1-1+l2+1-2*q;
                    B2 += a_l1ml2p[q]*I(self,nu,2*m);
                    //printf("a[%d] = %.15g (%d,%d)\n", q, (double)a_l1ml2m[q],l1-1,l2+1);
                }
                B2 *= -gaunt_a0(l1-1,l2+1,m);
                B2 *= (l1+1.)*(l1+m)*l2*(l2-m+1.);
            }

            B3 = 0;
            if((l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1pl2m; q++)
                {
                    nu = l1+1+l2-1-2*q;
                    B3 += a_l1pl2m[q]*I(self,nu,2*m);
                }
                B3 *= -gaunt_a0(l1+1,l2-1,m);
                B3 *= l1*(l1-m+1.)*(l2+1.)*(l2+m);
            }

            B4 = 0;
            for(q = 0; q <= qmax_l1pl2p; q++)
            {
                nu = l1+1+l2+1-2*q;
                B4 += a_l1pl2p[q]*I(self,nu,2*m);
            }
            B4 *= gaunt_a0(l1+1,l2+1,m);
            B4 *= l1*(l1-m+1.)*l2*(l2-m+1.);

            B = B1+B2+B3+B4;

            log_B  = logq(fabsq(B)) - logq(2*l1+1)-logq(2*l2+1);
            sign_B = copysignq(1,B);
        }

        /* C */
        {
            edouble C,C1,C2;

            C1 = 0;
            if((l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1l2m; q++)
                {
                    nu = l1+l2-1-2*q;
                    C1 += a_l1l2m[q]*I(self,nu,2*m);
                }
                C1 *= gaunt_a0(l1,l2-1,m);
                C1 *= (l2+1)*(l2+m);
            }

            C2 = 0;
            for(q = 0; q <= qmax_l1l2p; q++)
            {
                nu = l1+l2+1-2*q;
                C2 += a_l1l2p[q]*I(self,nu,2*m);
            }
            C2 *= -gaunt_a0(l1,l2+1,m);
            C2 *= l2*(l2-m+1);

            C = (C1+C2)/(2*l2+1);
            log_C  = logq(fabsq(C));
            sign_C = copysignq(1,C);
        }

        /* D */
        {
            edouble D,D1,D2;

            D1 = 0;
            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2; q++)
                {
                    nu = l1-1+l2-2*q;
                    D1 += a_l1ml2[q]*I(self,nu,2*m);
                }
                D1 *= gaunt_a0(l1-1,l2,m);
                D1 *= (l1+1)*(l1+m);
            }

            D2 = 0;
            for(q = 0; q <= qmax_l1pl2; q++)
            {
                nu = l1+1+l2-2*q;
                D2 += a_l1pl2[q]*I(self,nu,2*m);
            }
            D2 *= -gaunt_a0(l1+1,l2,m);
            D2 *= l1*(l1-m+1);

            D = (D1+D2)/(2*l1+1);
            log_D  = logq(fabsq(D));
            sign_D = copysignq(1,D);
        }


        cint->lnA_TM   = cint->lnA_TE = 2*log_m+lnLambda-tau+log_A;
        cint->signA_TM = -MPOW(l2+m)*sign_A; /* - because Lambda(l1,l2,m) is negative */
        cint->signA_TE = -cint->signA_TM;

        cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+log_B;
        cint->signB_TM = -MPOW(l2+m+1)*sign_B;
        cint->signB_TE = -cint->signB_TM;

        cint->lnC_TM   = cint->lnC_TE = log_m+lnLambda-tau+log_C;
        cint->signC_TM = -MPOW(l2+m)*sign_C;
        cint->signC_TE = -cint->signC_TM;

        cint->lnD_TM   = cint->lnD_TE = log_m+lnLambda-tau+log_D;
        cint->signD_TM = -MPOW(l2+m+1)*sign_D; // XXX ???
        cint->signD_TE = -cint->signD_TM;
    }
}

void cache_gaunt_init(gaunt_cache_t *cache, int lmax)
{
    int i;
    int N = lmax+2;
    size_t size = (N*N-N)/2+N;

    cache->N = N;
    cache->size = size;
    cache->cache = xmalloc(size*sizeof(edouble *));

    for(i = 0; i < size; i++)
        cache->cache[i] = NULL;
}

void cache_gaunt_free(gaunt_cache_t *cache)
{
    int i;
    size_t size = cache->size;

    for(i = 0; i < size; i++)
        if(cache->cache[i] != NULL)
            xfree(cache->cache[i]);

    xfree(cache->cache);
}

edouble *cache_gaunt_get(gaunt_cache_t *cache, int n, int nu, int m)
{
    int index;
    int N = cache->N;
    edouble *v;

    if(n <= nu)
        index = n*N - (n-1)*((n-1) + 1)/2 + nu - n;
    else
        index = nu*N - (nu-1)*((nu-1) + 1)/2 + n - nu;

    //printf("index=%d, N=%d\n", index, ((N*N-N)/2+N));
    v = cache->cache[index];
    if(v == NULL)
    {
        int len = gaunt_qmax(n,nu,m)+20; // XXX WTF?!?
        v = xmalloc(len*sizeof(edouble));
        gaunt(n, nu, m, v);
        cache->cache[index] = v;
    }

    return v;
}
