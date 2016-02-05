from __future__ import division, print_function
from mpmath import zeta
from math import log,pi,ceil,factorial

def H(q):
    """Return Harmonic number H(q)"""
    sum = 0
    for k in range(1,q+1):
        sum += 1/k
    return sum


def Lin(n,x,D=17):
    """Calculate polylogarithm of order n for 0 <= x <= 1 with D-decimal digit
    precision according to [1].

    References:
    [1] Crandall, Notes on fast polylogarithm computation, 2006
    """
    if x == 0:
        return 0
    elif x == 1:
        return zeta(n,1)

    L = int(ceil(D*log(10)/log(4)))

    if x < 0.25:
        sum = 0
        for k in range(1,L+1):
            sum += x**k/k**n
        return sum
    else:
        logx = log(x)
        sum = logx**(n-1)/factorial(n-1)*(H(n-1) - log(-logx))
        for m in range(L+1):
            if n-1 != m:
                sum += zeta(n-m,1)/factorial(m)*logx**m
        return sum


def Li2(x):
    """Implement the polylogarithm Li_s(x) for s=2 (also known as dilogarithm
    or Spencer function) for 0 <= x <= 1.

    This function uses a series with n^4 convergence, see [1]. In order to
    improve convergence, we make use of the identity:
        Li2(x)+Li2(1-x) = pi^2/6 - log(x)*log(1-x)

    References:
    [1] Robert Morris, The dilogarithm function of a real argument, Math. Comp. 33 (1979)
    """
    if x == 0:
        return 0
    elif x == 1:
        return pi**2/6
    elif x < 0.5:
        sum = 0
        for n in range(1,35):
            sum += x**n/(n*(n+1))**2
        return x/(x+1)*(3+sum)-2*(x-1)/(x+1)*log(1-x)
    else:
        return pi**2/6-log(x)*log(1-x)-Li2(1-x)
