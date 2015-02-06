#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sfunc.h"

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

    for(k = offset; k < len_p+offset; k++)
        //value += gamma(k+1)/tau**(k+1)*p[k-offset]
        value += expq(lgammaq(k+1)-(k+1)*logq(tau))*p[k-offset];

    return value;
}


static edouble I(int alpha, int beta, int nu, int m2, edouble tau)
{
    // exp(-z*tau) * (z^2+2z)^alpha * (1+z)^beta * Plm(nu, 2m, 1+z)
    int m = m2/2;
    edouble p1[m+alpha+beta+1], p2[nu+1-m2], p[-m+alpha+beta+1+nu];

    poly1(m+alpha, beta, p1);
    poly2(nu,m2,p2);
    polymult(p1, m+alpha+beta+1, p2, nu+1-m2, p);

    return MPOW(m)*polyintegrate(p, -m+alpha+beta+1+nu, m+alpha, tau);
}


edouble integralA(int l1, int l2, int m, edouble tau)
{
    const int qmax = GAUNT_QMAX(l1,l2,m,m);
    int q,nu;
    edouble prefactor, value = 0, a0 = GAUNT_a0(l1,l2,m,m);
    edouble a[qmax];

    if(m == 0)
        return 0;

    //prefactor = MPOW(l2)*pow_2(m)*expq(-tau);
    prefactor = 1;

    gaunt(l1, l2, m, m, a);
    for(q = 0; q <= qmax; q++)
    {
        nu = l1+l2-2*q;
        value += a[q]*I(-1,0,nu,2*m,tau);
    }

    return prefactor*a0*value;
}


void casimir_integrate_perf(casimir_integrals_t *cint, int l1, int l2, int m, double nT)
{
    edouble tau = 2*nT;
}

int main(int argc, char *argv[])
{
    int l1,l2,m;
    edouble tau;

    l1  = 4;
    l2  = 3;
    m   = 1;
    tau = 2;

    printf("IntegralA = %.14g\n", (double)integralA(l1,l2,m,tau));

    return 0;
}

/*

def integralB(l1,l2,m,tau):
    #prefactor = (-1)**(l2+1)*exp(-tau)
    prefactor = 1
    alpha = -1

    qmax,a0,a = gaunt(l1+1,l2+1,m,m)
    sum1 = 0
    for q in range(qmax+1):
        nu = l1+1+l2+1-2*q
        sum1 += a[q]*integral(alpha,0,nu,2*m,tau)
    sum1 *= a0*(l1-m+1)*(l2-m+1)

    qmax,a0,a = gaunt(l1+1,l2,m,m)
    sum2 = 0
    for q in range(qmax+1):
        nu = l1+1+l2-2*q
        sum2 += a[q]*integral(alpha,1,nu,2*m,tau)
    sum2 *= -a0*(l1-m+1)*(l2+1)

    qmax,a0,a = gaunt(l1,l2+1,m,m)
    sum3 = 0
    for q in range(qmax+1):
        nu = l1+1+l2-2*q
        sum3 += a[q]*integral(alpha,1,nu,2*m,tau)
    sum3 *= -a0*(l1+1)*(l2-m+1)

    qmax,a0,a = gaunt(l1,l2,m,m)
    sum4 = 0
    for q in range(qmax+1):
        nu = l1+l2-2*q
        sum4 += a[q]*integral(alpha,2,nu,2*m,tau)
    sum4 *= a0*(l1+1)*(l2+1)

    return prefactor*(sum1+sum2+sum3+sum4)


def integralC(l1,l2,m,tau):
    if m == 0:
        return 0

    #prefactor = (-1)**l2*m*exp(-tau)
    prefactor = 1

    qmax,a0,a = gaunt(l1,l2+1,m,m)
    sum1 = 0
    for q in range(qmax+1):
        nu = l1+l2+1-2*q
        sum1 += a[q]*integral(-1,0,nu,2*m,tau)
    sum1 *= a0*(l2-m+1)

    qmax,a0,a = gaunt(l1,l2,m,m)
    sum2 = 0
    for q in range(qmax+1):
        nu = l1+l2-2*q
        sum2 += a[q]*integral(-1,1,nu,2*m,tau)
    sum2 *= -a0*(l2+1)

    return prefactor*(sum1+sum2)


#print integralA(4,3,2,2)
#print integralB(4,3,2,2)
print integralB(4,3,1,2)
*/
