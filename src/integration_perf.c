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
    if(isinf(v))
    {
        edouble tau = self->tau;

        edouble p1[m], p2[nu+1-m2], p[-m+nu];

        poly1(m-1, p1);
        poly2(nu,m2,p2);
        polymult(p1, m, p2, nu+1-m2, p);

        v = polyintegrate(p, -m+nu, m-1, tau, NULL);
        //printf("v=%g\n", (double)v);
        self->cache_I[index] = v;

        TERMINATE(isinf(v), "I is inf, nu=%d, 2m=%d\n", nu, m2);
        TERMINATE(isnan(v), "I is nan, nu=%d, 2m=%d\n", nu, m2);
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
    int N = lmax+2;
    size_t elems = (N*N-N)/2+N;
    self->tau = 2*nT;

    /* init cache */
    self->nu_max = 2*lmax+4;
    self->m2_max = 1*lmax+2;
    size = self->nu_max*self->m2_max;
    self->cache_I = xmalloc(size*sizeof(edouble));

    for(i = 0; i < size; i++)
        self->cache_I[i] = -INFINITY;

    self->cache = NULL;
    self->signs = NULL;
    self->m     = -1;
    self->lmax  = lmax;

    self->N = N;
    self->elems = elems;
    self->cache = xmalloc(elems*sizeof(edouble *));
    self->signs = xmalloc(elems*sizeof(sign_t *));

    for(i = 0; i < elems; i++)
    {
        self->cache[i] = NULL;
        self->signs[i] = NULL;
    }
}

void casimir_integrate_perf_free(integration_perf_t *self)
{
    int i;
    for(i = 0; i < self->elems; i++)
    {
        if(self->cache[i] != NULL)
        {
            xfree(self->cache[i]);
            self->cache[i] = NULL;
        }
        if(self->signs[i] != NULL)
        {
            xfree(self->signs[i]);
            self->signs[i] = NULL;
        }
    }

    xfree(self->cache);
    self->cache = NULL;
    xfree(self->signs);
    self->signs = NULL;

    xfree(self->cache_I);
    self->cache_I = NULL;
}

/* TODO: m = 0, check signs! */
void casimir_integrate_perf(integration_perf_t *self, int l1, int l2, int m, casimir_integrals_t *cint)
{
    edouble tau = self->tau;
    edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);

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

        TERMINATE(isinf(cint->lnB_TM), "lnB is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        TERMINATE(isnan(cint->lnB_TM), "lnB is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
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

        edouble *sum = xmalloc((qmax_l1pl2p+1)*sizeof(edouble));

        sign_t *signs_l1l2;
        edouble *a_l1l2 = cache_gaunt_get(self, l1,l2,m, &signs_l1l2);

        sign_t *signs_l1pl2, *signs_l1ml2, *signs_l1l2p, *signs_l1l2m;
        edouble *a_l1pl2 = cache_gaunt_get(self, l1+1,l2,m, &signs_l1pl2);
        edouble *a_l1ml2 = cache_gaunt_get(self, l1-1,l2,m, &signs_l1ml2);
        edouble *a_l1l2p = cache_gaunt_get(self, l1,l2+1,m, &signs_l1l2p);
        edouble *a_l1l2m = cache_gaunt_get(self, l1,l2-1,m, &signs_l1l2m);

        sign_t *signs_l1pl2p, *signs_l1pl2m, *signs_l1ml2p, *signs_l1ml2m;
        edouble *a_l1pl2p = cache_gaunt_get(self, l1+1,l2+1,m, &signs_l1pl2p);
        edouble *a_l1pl2m = cache_gaunt_get(self, l1+1,l2-1,m, &signs_l1pl2m);
        edouble *a_l1ml2p = cache_gaunt_get(self, l1-1,l2+1,m, &signs_l1ml2p);
        edouble *a_l1ml2m = cache_gaunt_get(self, l1-1,l2-1,m, &signs_l1ml2m);


        /* A */
        {
            for(q = 0; q <= qmax_l1l2; q++)
            {
                nu = l1+l2-2*q;
                sum[q] = a_l1l2[q] + I(self,nu,2*m);
            }

            log_A = gaunt_log_a0(l1,l2,m)+logadd_ms(sum, signs_l1l2, qmax_l1l2+1, &sign_A);
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
                    sum[q] = a_l1ml2m[q]+I(self,nu,2*m);
                }
                log_Bn[0]  = gaunt_log_a0(l1-1,l2-1,m) + logadd_ms(sum, signs_l1ml2m, qmax_l1ml2m+1, &signs[0]);
                log_Bn[0] += logq(l1+1)+logq(l1+m)+logq(l2+1)+logq(l2+m);
            }

            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2p; q++)
                {
                    nu = l1-1+l2+1-2*q;
                    sum[q] = a_l1ml2p[q] + I(self,nu,2*m);
                }
                log_Bn[1]  = gaunt_log_a0(l1-1,l2+1,m) + logadd_ms(sum, signs_l1ml2p, qmax_l1ml2p+1, &signs[1]);
                log_Bn[1] += logq(l1+1)+logq(l1+m)+logq(l2)+logq(l2-m+1);
                signs[1] *= -1;
            }

            if((l2-1) >= m)
            {
                for(q = 0; q <= qmax_l1pl2m; q++)
                {
                    nu = l1+1+l2-1-2*q;
                    sum[q] = a_l1pl2m[q] + I(self,nu,2*m);
                }
                log_Bn[2]  = gaunt_log_a0(l1+1,l2-1,m) + logadd_ms(sum, signs_l1pl2m, qmax_l1pl2m+1, &signs[2]);
                log_Bn[2] += logq(l1)+logq(l1-m+1)+logq(l2+1)+logq(l2+m);
                signs[2] *= -1;
            }

            for(q = 0; q <= qmax_l1pl2p; q++)
            {
                nu = l1+1+l2+1-2*q;
                sum[q] = a_l1pl2p[q] + I(self,nu,2*m);
            }
            log_Bn[3]  = gaunt_log_a0(l1+1,l2+1,m) + logadd_ms(sum, signs_l1pl2p, qmax_l1pl2p+1, &signs[3]);
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
                    sum[q] = a_l1l2m[q] + I(self,nu,2*m);
                }
                log_Cn[0]  = gaunt_log_a0(l1,l2-1,m) + logadd_ms(sum, signs_l1l2m, qmax_l1l2m+1, &signs[0]);
                log_Cn[0] += logq(l2+1)+logq(l2+m);
            }

            for(q = 0; q <= qmax_l1l2p; q++)
            {
                nu = l1+l2+1-2*q;
                sum[q] = a_l1l2p[q] + I(self,nu,2*m);
            }
            log_Cn[1] = gaunt_log_a0(l1,l2+1,m) + logadd_ms(sum, signs_l1l2p, qmax_l1l2p+1, &signs[1]);
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
                    sum[q] = a_l1ml2[q] + I(self,nu,2*m);
                }
                log_Dn[0]  = gaunt_log_a0(l1-1,l2,m) + logadd_ms(sum, signs_l1ml2, qmax_l1ml2+1, &signs[0]);
                log_Dn[0] += logq(l1+1)+logq(l1+m);
            }

            for(q = 0; q <= qmax_l1pl2; q++)
            {
                nu = l1+1+l2-2*q;
                sum[q] = a_l1pl2[q] + I(self,nu,2*m);
            }
            log_Dn[1]  = gaunt_log_a0(l1+1,l2,m) + logadd_ms(sum, signs_l1pl2, qmax_l1pl2+1, &signs[1]);
            log_Dn[1] += logq(l1)+logq(l1-m+1);
            signs[1] *= -1;

            log_D = logadd_ms(log_Dn, signs, 2, &sign_D) - logq(2*l1+1);
        }

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

        TERMINATE(isinf(cint->lnA_TM), "lnA is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        TERMINATE(isnan(cint->lnA_TM), "lnA is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

        TERMINATE(isinf(cint->lnB_TM), "lnB is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        TERMINATE(isnan(cint->lnB_TM), "lnB is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

        TERMINATE(isinf(cint->lnC_TM), "lnC is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        TERMINATE(isnan(cint->lnC_TM), "lnC is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);

        TERMINATE(isinf(cint->lnD_TM), "lnD is inf, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
        TERMINATE(isnan(cint->lnD_TM), "lnD is nan, l1=%d,l2=%d,m=%d,tau=%g", l1,l2,m,(double)tau);
    }
}

edouble *cache_gaunt_get(integration_perf_t *self, int n, int nu, int m, sign_t **signs)
{
    int index;
    const int N = self->N;
    edouble *v;

    if(self->m != m)
    {
        int i;
        for(i = 0; i < self->elems; i++)
        {
            if(self->cache[i] != NULL)
            {
                xfree(self->cache[i]);
                self->cache[i] = NULL;
            }
            if(self->signs[i] != NULL)
            {
                xfree(self->signs[i]);
                self->signs[i] = NULL;
            }
        }
        self->m = m;
    }

    #define INDEX1(n,nu) ((n)*(N) - ((n)-1)*(((n)-1) + 1)/2 + (nu) - (n))
    #define INDEX2(n,nu) ((nu)*(N) - ((nu)-1)*(((nu)-1) + 1)/2 + (n) - (nu))

    if(n <= nu)
        index = INDEX1(n,nu);
    else
        index = INDEX2(n,nu);

    v = self->cache[index];
    if(v == NULL)
    {
        int i;
        int elems = MAX(0,gaunt_qmax(n,nu,m))+1;

        /* this could be improved */
        if(n > 3)
        {
            int q;
            for(q = 0; q < n-3; q++)
            {
                int index_new = INDEX2(n-3,q);
                if(self->cache[index_new] != NULL)
                {
                    xfree(self->cache[index_new]);
                    self->cache[index_new] = NULL;
                }
                if(self->cache[index_new] != NULL)
                {
                    xfree(self->cache[index_new]);
                    self->cache[index_new] = NULL;
                }
            }
        }

        self->cache[index] = xmalloc(elems*sizeof(edouble));
        self->signs[index] = xmalloc(elems*sizeof(sign_t));
        gaunt(n, nu, m, self->cache[index]);
        for(i = 0; i < elems; i++)
        {
            self->signs[index][i] = copysignq(1, self->cache[index][i]);
            self->cache[index][i] = logq(fabsq(self->cache[index][i]));
        }

        *signs = self->signs[index];
        return self->cache[index];
    }
    else
    {
        *signs = self->signs[index];
        return v;
    }
}
