#include <assert.h>
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


static edouble polyintegrate(edouble p[], int len_p, int offset, edouble tau, sign_t *sign)
{
    int k;
    edouble value = 0;
    edouble max = lgammaq(len_p+offset);

    for(k = offset; k < len_p+offset; k++)
        value += expq(lgammaq(k+1)-(k+1)*logq(tau)-max)*p[k-offset];

    if(sign != NULL)
        *sign = copysignq(1,value);

    return logq(fabsq(value)) + max;
}


/* evaluete integral I_nu^2m(tau) = (-1)^m * exp(-z*tau)/ (z^2+2z) * Plm(nu, 2m, 1+z) */
static inline edouble I(integration_perf_t *self, int nu, int m2)
{
    const int m = m2/2;
    const int index = m*self->nu_max + nu;
    edouble v;

    v = self->cache_I[index];
    if(isinfq(v))
    {
        edouble tau = self->tau;

        edouble p1[m], p2[nu+1-m2], p[-m+nu];

        poly1(m-1, p1);
        poly2(nu,m2,p2);
        polymult(p1, m, p2, nu+1-m2, p);

        v = polyintegrate(p, -m+nu, m-1, tau, NULL);
        //printf("v=%g\n", (double)v);
        self->cache_I[index] = v;

        #ifndef NDEBUG
            if(isinfq(v))
                TERMINATE("I is inf, nu=%d, 2m=%d\n", nu, m2);
            if(isnanq(v))
                TERMINATE("I is nan, nu=%d, 2m=%d\n", nu, m2);
        #endif
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
        self->cache_I[i] = -INFINITY;

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
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint)
{
    edouble tau = self->tau;
    edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);
    gaunt_cache_t *cache = self->cache_gaunt;

    if(m == 0)
    {
        sign_t sign;
        edouble value;
        edouble pm[] = {0,2,1}; /* z²+2z */
        edouble pdpl1m[l1-m+2];
        edouble pdpl2m[l2-m+2];
        edouble pdpl1mpdpl2m[l1+l2-1];
        edouble pmpdpl1mpdpl2m[l1+l2+1];

        polydplm(pdpl1m,pdpl2m,l1,l2,m);

        polymult(pdpl1m, l1, pdpl2m, l2, pdpl1mpdpl2m);
        polymult(pm, 3, pdpl1mpdpl2m, l1+l2-1, pmpdpl1mpdpl2m);

        value = polyintegrate(pmpdpl1mpdpl2m, l1+l2+1, 0, tau, &sign);

        cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+value;
        cint->signB_TM = -MPOW(l2+1)*sign;
        cint->signB_TE = -cint->signB_TM;

        cint->lnA_TM = cint->lnA_TE = -INFINITY;
        cint->lnC_TM = cint->lnC_TE = -INFINITY;
        cint->lnD_TM = cint->lnD_TE = -INFINITY;
        cint->signA_TM = cint->signA_TE = +1;
        cint->signC_TM = cint->signC_TE = +1;
        cint->signD_TM = cint->signD_TE = +1;

        #ifndef NDEBUG
            if(isinfq(cint->lnB_TM))
                TERMINATE("lnB is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
            if(isnanq(cint->lnB_TM))
                TERMINATE("lnB is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        #endif
    }
    else
    {
        int nu,q;
        sign_t sign_A, sign_B, sign_C, sign_D;
        const edouble log_m = logq(m);
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

        edouble *sum;
        sign_t *signs_sum;

        edouble *a_l1l2 = cache_gaunt_get(cache, l1,l2,m);

        edouble *a_l1pl2 = cache_gaunt_get(cache, l1+1,l2,m);
        edouble *a_l1ml2 = cache_gaunt_get(cache, l1-1,l2,m);
        edouble *a_l1l2p = cache_gaunt_get(cache, l1,l2+1,m);
        edouble *a_l1l2m = cache_gaunt_get(cache, l1,l2-1,m);

        edouble *a_l1pl2p = cache_gaunt_get(cache, l1+1,l2+1,m);
        edouble *a_l1pl2m = cache_gaunt_get(cache, l1+1,l2-1,m);
        edouble *a_l1ml2p = cache_gaunt_get(cache, l1-1,l2+1,m);
        edouble *a_l1ml2m = cache_gaunt_get(cache, l1-1,l2-1,m);

        sum       = xmalloc((qmax_l1pl2p+1)*sizeof(edouble));
        signs_sum = xmalloc((qmax_l1pl2p+1)*sizeof(sign_t));

        /* A */
        {
            for(q = 0; q <= qmax_l1l2; q++)
            {
                nu = l1+l2-2*q;
                sum[q]       = logq(fabsq(a_l1l2[q])) + I(self,nu,2*m);
                signs_sum[q] = copysignq(1,a_l1l2[q]);
            }

            log_A = gaunt_log_a0(l1,l2,m)+logadd_ms(sum, signs_sum, qmax_l1l2+1, &sign_A);
        }

        /* B */
        {
            sign_t signs[4] = {1,1,1,1};
            edouble log_Bn[4] = { -INFINITY, -INFINITY, -INFINITY, -INFINITY };

            if((l1-1) >= m && (l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2m; q++)
                {
                    nu = l1-1+l2-1-2*q;
                    sum[q] = logq(fabsq(a_l1ml2m[q]))+I(self,nu,2*m);
                    signs_sum[q] = copysignq(1, a_l1ml2m[q]);
                }
                log_Bn[0]  = gaunt_log_a0(l1-1,l2-1,m) + logadd_ms(sum, signs_sum, qmax_l1ml2m+1, &signs[0]);
                log_Bn[0] += logq(l1+1)+logq(l1+m)+logq(l2+1)+logq(l2+m);
            }

            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2p; q++)
                {
                    nu = l1-1+l2+1-2*q;
                    sum[q] = logq(fabsq(a_l1ml2p[q]))+ + I(self,nu,2*m);
                    signs_sum[q] = copysignq(1, a_l1ml2p[q]);
                }
                log_Bn[1]  = gaunt_log_a0(l1-1,l2+1,m) + logadd_ms(sum, signs_sum, qmax_l1ml2p+1, &signs[1]);
                log_Bn[1] += logq(l1+1)+logq(l1+m)+logq(l2)+logq(l2-m+1);
                signs[1] *= -1;
            }

            if((l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1pl2m; q++)
                {
                    nu = l1+1+l2-1-2*q;
                    sum[q] = logq(fabsq(a_l1pl2m[q])) + I(self,nu,2*m);
                    signs_sum[q] = copysignq(1,a_l1pl2m[q]);
                }
                log_Bn[2]  = gaunt_log_a0(l1+1,l2-1,m) + logadd_ms(sum, signs_sum, qmax_l1pl2m+1, &signs[2]);
                log_Bn[2] += logq(l1)+logq(l1-m+1)+logq(l2+1)+logq(l2+m);
                signs[2] *= -1;
            }

            for(q = 0; q <= qmax_l1pl2p; q++)
            {
                nu = l1+1+l2+1-2*q;
                sum[q] = logq(fabsq(a_l1pl2p[q])) + I(self,nu,2*m);
                signs_sum[q] = copysignq(1,a_l1pl2p[q]);
            }
            log_Bn[3]  = gaunt_log_a0(l1+1,l2+1,m) + logadd_ms(sum, signs_sum, qmax_l1pl2p+1, &signs[3]);
            log_Bn[3] += logq(l1)+logq(l1-m+1)+logq(l2)+logq(l2-m+1);

            log_B = logadd_ms(log_Bn, signs, 4, &sign_B) - logq(2*l1+1) - logq(2*l2+1);
        }

        /* C */
        {
            sign_t signs[2] = {1,1};
            edouble log_Cn[2] = { -INFINITY, -INFINITY };

            if((l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1l2m; q++)
                {
                    nu = l1+l2-1-2*q;
                    sum[q]       = logq(fabsq(a_l1l2m[q])) + I(self,nu,2*m);
                    signs_sum[q] = copysignq(1, a_l1l2m[q]);
                }
                log_Cn[0]  = gaunt_log_a0(l1,l2-1,m) + logadd_ms(sum, signs_sum, qmax_l1l2m+1, &signs[0]);
                log_Cn[0] += logq(l2+1)+logq(l2+m);
            }

            for(q = 0; q <= qmax_l1l2p; q++)
            {
                nu = l1+l2+1-2*q;
                sum[q] = logq(fabsq(a_l1l2p[q])) + I(self,nu,2*m);
                signs_sum[q] = copysignq(1, a_l1l2p[q]);
            }
            log_Cn[1] = gaunt_log_a0(l1,l2+1,m) + logadd_ms(sum, signs_sum, qmax_l1l2p+1, &signs[1]);
            log_Cn[1] += logq(l2)+logq(l2-m+1);
            signs[1] *= -1;

            log_C = logadd_ms(log_Cn, signs, 2, &sign_C) - logq(2*l2+1);
        }

        /* D */
        {
            sign_t signs[2] = {1,1};
            edouble log_Dn[2] = { -INFINITY, -INFINITY };

            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2; q++)
                {
                    nu = l1-1+l2-2*q;
                    sum[q] = logq(fabsq(a_l1ml2[q])) + I(self,nu,2*m);
                    signs_sum[q] = copysignq(1,a_l1ml2[q]);
                }
                log_Dn[0]  = gaunt_log_a0(l1-1,l2,m) + logadd_ms(sum, signs_sum, qmax_l1ml2+1, &signs[0]);
                log_Dn[0] += logq(l1+1)+logq(l1+m);
            }

            for(q = 0; q <= qmax_l1pl2; q++)
            {
                nu = l1+1+l2-2*q;
                sum[q] = logq(fabsq(a_l1pl2[q])) + I(self,nu,2*m);
                signs_sum[q] = copysignq(1, a_l1pl2[q]);
            }
            log_Dn[1]  = gaunt_log_a0(l1+1,l2,m) + logadd_ms(sum, signs_sum, qmax_l1pl2+1, &signs[1]);
            log_Dn[1] += logq(l1)+logq(l1-m+1);
            signs[1] *= -1;

            log_D = logadd_ms(log_Dn, signs, 2, &sign_D) - logq(2*l1+1);
        }

        xfree(signs_sum);
        xfree(sum);


        cint->lnA_TM   = cint->lnA_TE = 2*log_m+lnLambda-tau+log_A;
        cint->signA_TM = -MPOW(l2)*sign_A; /* - because Lambda(l1,l2,m) is negative */
        cint->signA_TE = -cint->signA_TM;

        cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+log_B;
        cint->signB_TM = -MPOW(l2+1)*sign_B;
        cint->signB_TE = -cint->signB_TM;

        cint->lnC_TM   = cint->lnC_TE = log_m+lnLambda-tau+log_C;
        cint->signC_TM = -MPOW(l2)*sign_C;
        cint->signC_TE = -cint->signC_TM;

        cint->lnD_TM   = cint->lnD_TE = log_m+lnLambda-tau+log_D;
        cint->signD_TM = -MPOW(l2+1)*sign_D;
        cint->signD_TE = -cint->signD_TM;

        #ifndef NDEBUG
            if(isinfq(cint->lnA_TM))
                TERMINATE("lnA is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
            if(isnanq(cint->lnA_TM))
                TERMINATE("lnA is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

            if(isinfq(cint->lnB_TM))
                TERMINATE("lnB is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
            if(isnanq(cint->lnB_TM))
                TERMINATE("lnB is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

            if(isinfq(cint->lnC_TM))
                TERMINATE("lnC is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
            if(isnanq(cint->lnC_TM))
                TERMINATE("lnC is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

            if(isinfq(cint->lnD_TM))
                TERMINATE("lnD is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
            if(isnanq(cint->lnD_TM))
                TERMINATE("lnD is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        #endif
    }
}

void cache_gaunt_init(gaunt_cache_t *cache, int lmax)
{
    int i;
    int N = lmax+2;
    size_t elems = (N*N-N)/2+N;

    cache->N = N;
    cache->elems = elems;
    cache->cache = xmalloc(elems*sizeof(edouble *));
    cache->signs = xmalloc(elems*sizeof(sign_t *));

    for(i = 0; i < elems; i++)
    {
        cache->cache[i] = NULL;
        cache->signs[i] = NULL;
    }
}

void cache_gaunt_free(gaunt_cache_t *cache)
{
    int i;
    size_t elems = cache->elems;

    for(i = 0; i < elems; i++)
    {
        if(cache->cache[i] != NULL)
        {
            xfree(cache->cache[i]);
            cache->cache[i] = NULL;
        }
        if(cache->signs[i] != NULL)
        {
            xfree(cache->signs[i]);
            cache->signs[i] = NULL;
        }
    }

    if(cache->cache != NULL)
    {
        xfree(cache->cache);
        cache->cache = NULL;
    }
    if(cache->signs != NULL)
    {
        xfree(cache->signs);
        cache->signs = NULL;
    }
}

edouble *cache_gaunt_get(gaunt_cache_t *cache, int n, int nu, int m)
{
    int index;
    const int N = cache->N;
    edouble *v;

    #define INDEX1(n,nu) ((n)*(N) - ((n)-1)*(((n)-1) + 1)/2 + (nu) - (n))
    #define INDEX2(n,nu) ((nu)*(N) - ((nu)-1)*(((nu)-1) + 1)/2 + (n) - (nu))

    if(n <= nu)
        index = INDEX1(n,nu);
    else
        index = INDEX2(n,nu);

    v = cache->cache[index];
    if(v == NULL)
    {
        size_t elems = MAX(0,gaunt_qmax(n,nu,m))+1;

        /* this could be improved */
        if(n > 3)
        {
            int q;
            for(q = 0; q < n-3; q++)
            {
                int index_new = INDEX2(n-3,q);
                if(cache->cache[index_new] != NULL)
                {
                    xfree(cache->cache[index_new]);
                    cache->cache[index_new] = NULL;
                }
                if(cache->cache[index_new] != NULL)
                {
                    xfree(cache->cache[index_new]);
                    cache->cache[index_new] = NULL;
                }
            }
        }

        cache->cache[index] = xmalloc(elems*sizeof(edouble));
        gaunt(n, nu, m, cache->cache[index]);
        return cache->cache[index];
    }
    else
        return v;
}
