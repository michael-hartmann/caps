/**
 * @file   integration_perf.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   August, 2015
 * @brief  Perform integration for perfect reflectors
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "edouble.h"
#include "sfunc.h"
#include "libcasimir.h"
#include "integration_perf.h"

#include "utils.h"


static edouble polyintegrate(edouble p[], const int len_p, const int offset, const edouble tau)
{
    int k;
    const edouble log_tau = loge(tau);
    edouble list[len_p];

    for(k = offset; k < len_p+offset; k++)
        list[k-offset] = lgammae(k+1)-(k+1)*log_tau+p[k-offset];

    return logadd_m(list, len_p);
}


static void log_polymult(edouble p1[], const int len_p1, edouble p2[], const int len_p2, edouble p[])
{
    int k, len = len_p1+len_p2-1;
    edouble temp[len];

    for(k = 0; k < len; k++)
    {
        int i, elems = 0;

        for(i = 0; i < MIN(k+1, len_p1); i++)
        {
            int j = k-i;
            if(i < len_p1 && j < len_p2)
                temp[elems++] = p1[i]+p2[j];
        }

        if(elems)
            p[k] = logadd_m(temp, elems);
        else
            p[k] = -INFINITY;
    }
}


/* evaluate integral I_nu^2m(tau) = (-1)^m * exp(-z*tau)/ (z^2+2z) * Plm(nu, 2m, 1+z) */
edouble _I(int nu, int m2, edouble tau)
{
    const int m = m2/2;
    edouble v;
    int k;

    edouble p1[m];       /* polynom (z+2)^(m-1) */
    edouble p2[nu+1-m2]; /* polynom d^(2m)/dz^(2m) P_(nu)(1+z) */
    edouble p[-m+nu];    /* polynom p1*p2 */

    /* Every monom of both polynoms is positive. So we can save the
     * logarithms of the coefficients. */

    for(k = 0; k <= m-1; k++)
        p1[k] = lgammae(m)-lgammae(k+1)-lgammae(m-k)+(m-1-k)*LOG2;

    for(k = m2; k <= nu; k++)
        p2[k-m2] = lgammae(k+nu+1)-lgammae(k+1)-lgammae(k-m2+1)-lgammae(-k+nu+1)-k*LOG2;

    log_polymult(p1, m, p2, nu+1-m2, p); /* len: nu-m */

    v = polyintegrate(p, -m+nu, m-1, tau);
    TERMINATE(!isfinite(v), "I=%Lg, nu=%d, 2m=%d\n", v, nu, m2);
    return v;
}

edouble I(integration_perf_t *self, int nu, int m2)
{
    const int m = m2/2;
    const int index = m*self->nu_max + nu;
    edouble v = self->cache_I[index];

    if(isinf(v))
        v = self->cache_I[index] = _I(nu, m2, self->tau);

    return v;
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
    int nu,q;
    sign_t sign_A, sign_B, sign_C, sign_D;
    edouble tau = self->tau;
    edouble lnLambda = casimir_lnLambda(l1, l2, m, NULL);
    edouble log_A,log_B,log_C,log_D, *sum = NULL;
    sign_t *signs_l1pl2p, *signs_l1pl2m, *signs_l1ml2p, *signs_l1ml2m;

    if(m == 0)
    {
        const int qmax_l1pl2p = gaunt_qmax(l1+1,l2+1,2);
        const int qmax_l1pl2m = gaunt_qmax(l1+1,l2-1,2);
        const int qmax_l1ml2p = gaunt_qmax(l1-1,l2+1,2);
        const int qmax_l1ml2m = gaunt_qmax(l1-1,l2-1,2);

        sum = xmalloc((qmax_l1pl2p+1)*sizeof(edouble));

        edouble *a_l1pl2p = cache_gaunt_get(self, l1+1,l2+1,2, &signs_l1pl2p);
        edouble *a_l1pl2m = cache_gaunt_get(self, l1+1,l2-1,2, &signs_l1pl2m);
        edouble *a_l1ml2p = cache_gaunt_get(self, l1-1,l2+1,2, &signs_l1ml2p);
        edouble *a_l1ml2m = cache_gaunt_get(self, l1-1,l2-1,2, &signs_l1ml2m);
        sign_t signs[4] = {1,1,1,1};
        edouble log_Bn[4] = { -INFINITY, -INFINITY, -INFINITY, -INFINITY };

        if((l1-1) >= 2 && (l2-1) >= 2)
        {
            for(q = 0; q <= qmax_l1ml2m; q++)
            {
                nu = l1-1+l2-1-2*q;
                sum[q] = a_l1ml2m[q]+I(self,nu,4);
            }
            log_Bn[0] = gaunt_log_a0(l1-1,l2-1,2) + logadd_ms(sum, signs_l1ml2m, qmax_l1ml2m+1, &signs[0]);
        }

        if((l1-1) >= 2)
        {
            for(q = 0; q <= qmax_l1ml2p; q++)
            {
                nu = l1-1+l2+1-2*q;
                sum[q] = a_l1ml2p[q] + I(self,nu,4);
            }
            log_Bn[1] = gaunt_log_a0(l1-1,l2+1,2) + logadd_ms(sum, signs_l1ml2p, qmax_l1ml2p+1, &signs[1]);
            signs[1] *= -1;
        }

        if((l2-1) >= 2)
        {
            for(q = 0; q <= qmax_l1pl2m; q++)
            {
                nu = l1+1+l2-1-2*q;
                sum[q] = a_l1pl2m[q] + I(self,nu,4);
            }
            log_Bn[2] = gaunt_log_a0(l1+1,l2-1,2) + logadd_ms(sum, signs_l1pl2m, qmax_l1pl2m+1, &signs[2]);
            signs[2] *= -1;
        }

        for(q = 0; q <= qmax_l1pl2p; q++)
        {
            nu = l1+1+l2+1-2*q;
            sum[q] = a_l1pl2p[q] + I(self,nu,4);
        }
        log_Bn[3] = gaunt_log_a0(l1+1,l2+1,2) + logadd_ms(sum, signs_l1pl2p, qmax_l1pl2p+1, &signs[3]);

        log_B = logadd_ms(log_Bn, signs, 4, &sign_B) - loge(2*l1+1) - loge(2*l2+1);

        cint->lnB_TM   = cint->lnB_TE = lnLambda-tau+log_B;
        cint->signB_TM = -MPOW(l2+1)*sign_B;
        cint->signB_TE = -cint->signB_TM;
        TERMINATE(!isfinite(cint->lnB_TM), "lnB=%Lg, l1=%d,l2=%d,m=%d,tau=%Lg", cint->lnB_TM,l1,l2,m,tau);

        cint->lnA_TE = cint->lnA_TM   = -INFINITY;
        cint->signA_TM = cint->signA_TE = 1;

        cint->lnC_TE = cint->lnC_TM   = -INFINITY;
        cint->signC_TM = cint->signC_TE = 1;

        cint->lnD_TE = cint->lnD_TM   = -INFINITY;
        cint->signD_TM = cint->signD_TE = 1;
    }
    else
    {
        const edouble log_m = loge(m);

        const int qmax_l1l2   = gaunt_qmax(l1,  l2,  m);

        const int qmax_l1pl2  = gaunt_qmax(l1+1,l2,  m);
        const int qmax_l1ml2  = gaunt_qmax(l1-1,l2,  m);
        const int qmax_l1l2p  = gaunt_qmax(l1,  l2+1,m);
        const int qmax_l1l2m  = gaunt_qmax(l1,  l2-1,m);

        const int qmax_l1pl2p = gaunt_qmax(l1+1,l2+1,m);
        const int qmax_l1pl2m = gaunt_qmax(l1+1,l2-1,m);
        const int qmax_l1ml2p = gaunt_qmax(l1-1,l2+1,m);
        const int qmax_l1ml2m = gaunt_qmax(l1-1,l2-1,m);

        sum = xmalloc((qmax_l1pl2p+1)*sizeof(edouble));

        sign_t *signs_l1l2;
        edouble *a_l1l2 = cache_gaunt_get(self, l1,l2,m, &signs_l1l2);

        sign_t *signs_l1pl2, *signs_l1ml2, *signs_l1l2p, *signs_l1l2m;
        edouble *a_l1pl2 = cache_gaunt_get(self, l1+1,l2,m, &signs_l1pl2);
        edouble *a_l1ml2 = cache_gaunt_get(self, l1-1,l2,m, &signs_l1ml2);
        edouble *a_l1l2p = cache_gaunt_get(self, l1,l2+1,m, &signs_l1l2p);
        edouble *a_l1l2m = cache_gaunt_get(self, l1,l2-1,m, &signs_l1l2m);

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
                log_Bn[0] += loge(l1+1)+loge(l1+m)+loge(l2+1)+loge(l2+m);
            }

            if((l1-1) >= m)
            {
                for(q = 0; q <= qmax_l1ml2p; q++)
                {
                    nu = l1-1+l2+1-2*q;
                    sum[q] = a_l1ml2p[q] + I(self,nu,2*m);
                }
                log_Bn[1]  = gaunt_log_a0(l1-1,l2+1,m) + logadd_ms(sum, signs_l1ml2p, qmax_l1ml2p+1, &signs[1]);
                log_Bn[1] += loge(l1+1)+loge(l1+m)+loge(l2)+loge(l2-m+1);
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
                log_Bn[2] += loge(l1)+loge(l1-m+1)+loge(l2+1)+loge(l2+m);
                signs[2] *= -1;
            }

            for(q = 0; q <= qmax_l1pl2p; q++)
            {
                nu = l1+1+l2+1-2*q;
                sum[q] = a_l1pl2p[q] + I(self,nu,2*m);
            }
            log_Bn[3]  = gaunt_log_a0(l1+1,l2+1,m) + logadd_ms(sum, signs_l1pl2p, qmax_l1pl2p+1, &signs[3]);
            log_Bn[3] += loge(l1)+loge(l1-m+1)+loge(l2)+loge(l2-m+1);

            log_B = logadd_ms(log_Bn, signs, 4, &sign_B) - loge(2*l1+1) - loge(2*l2+1);
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
                log_Cn[0] += loge(l2+1)+loge(l2+m);
            }

            for(q = 0; q <= qmax_l1l2p; q++)
            {
                nu = l1+l2+1-2*q;
                sum[q] = a_l1l2p[q] + I(self,nu,2*m);
            }
            log_Cn[1] = gaunt_log_a0(l1,l2+1,m) + logadd_ms(sum, signs_l1l2p, qmax_l1l2p+1, &signs[1]);
            log_Cn[1] += loge(l2)+loge(l2-m+1);
            signs[1] *= -1;

            log_C = logadd_ms(log_Cn, signs, 2, &sign_C) - loge(2*l2+1);
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
                log_Dn[0] += loge(l1+1)+loge(l1+m);
            }

            for(q = 0; q <= qmax_l1pl2; q++)
            {
                nu = l1+1+l2-2*q;
                sum[q] = a_l1pl2[q] + I(self,nu,2*m);
            }
            log_Dn[1]  = gaunt_log_a0(l1+1,l2,m) + logadd_ms(sum, signs_l1pl2, qmax_l1pl2+1, &signs[1]);
            log_Dn[1] += loge(l1)+loge(l1-m+1);
            signs[1] *= -1;

            log_D = logadd_ms(log_Dn, signs, 2, &sign_D) - loge(2*l1+1);
        }

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

        TERMINATE(!isfinite(cint->lnA_TM), "lnA=%Lg, l1=%d,l2=%d,m=%d,tau=%Lg", cint->lnA_TM,l1,l2,m,tau);
        TERMINATE(!isfinite(cint->lnB_TM), "lnB=%Lg, l1=%d,l2=%d,m=%d,tau=%Lg", cint->lnB_TM,l1,l2,m,tau);
        TERMINATE(!isfinite(cint->lnC_TM), "lnC=%Lg, l1=%d,l2=%d,m=%d,tau=%Lg", cint->lnC_TM,l1,l2,m,tau);
        TERMINATE(!isfinite(cint->lnD_TM), "lnD=%Lg, l1=%d,l2=%d,m=%d,tau=%Lg", cint->lnD_TM,l1,l2,m,tau);
    }

    xfree(sum);
}

/* We _do_ need a cache for the Gaunt coefficients or we just calculate and
 * recalculate the coefficients. However, this cache is crap. It works, but it
 * could be definitely improved. */
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
        size_t elems = MAX(0,1+gaunt_qmax(n,nu,m));

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
            self->signs[index][i] = copysigne(1, self->cache[index][i]);
            self->cache[index][i] = loge(fabse(self->cache[index][i]));
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
