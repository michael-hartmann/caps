import numpy as np
from math import sqrt, log1p,log10
from scipy.special import zeta
from scipy.optimize import bisect

zeta3 = zeta(3)


def E_HT_exact(LbyR, eps=1e-20):
    """Compute the exact Casimir interaction energy in the high-temperature
    limit for Drude metals. The value returned is in values of kb*T.

    Parameters:
        LbyR: L/R
        eps:  truncation of sum

    Return value:
        E/(kb*T)
    """
    x = LbyR
    Z = (1+x)*(1-sqrt((x**2+2*x)/(x**2+2*x+1)))

    sum1, sum2 = 0, 0
    l = 1
    while True:
        a = Z**(2*l) # Z^(2l)
        b = a*Z      # Z^(2l+1)
        v1 = (2*l+1)*log1p(-b)
        v2 = b*(1-a)/(1-b)

        sum1 += v1
        sum2 += v2

        if abs(v1/sum1) < eps and abs(v2/sum2) < eps:
            break

        l += 1

    return (sum1 + log1p(-(1-Z**2)*sum2))/2


def F_HT_exact(LbyR, eps=1e-20):
    """Compute the exact Casimir force in the high-temperature limit for Drude
    metals. The value returned is in values of kb*T/R.

    Parameters:
        LbyR: L/R
        eps:  truncation of sum

    Return value:
        F/(kb*T/R)
    """

    x = LbyR
    Z = (1+x)*(1-sqrt((x**2+2*x)/(x**2+2*x+1)))
    dZ = (sqrt(x*(2+x))-x-1)/sqrt(x*(2+x)) # dZ/dx

    sum1,sum2,sum3 = 0,0,0
    l = 1
    while True:
        a = Z**(2*l) # Z^(2l)
        b = a*Z      # Z^(2l+1)

        v1 = -(2*l+1)**2*a/(1-b)
        v2 = b*(1-a)/(1-b)
        v3 = a*( (a-1)*(1-3*Z**2+2*Z**(2*l+3)) + 2*l*(Z**2-1)*(Z**(4*l+1)-2*Z**(2*l)+1) )/(b-1)**2

        sum1 += v1
        sum2 += v2
        sum3 += v3

        if abs(v1/sum1) < eps and abs(v2/sum2) < eps and abs(v3/sum3):
            break

        l += 1

    return -(sum1 + sum3/(1-(1-Z**2)*sum2))*dZ/2


def F_HT_PFA(LbyR):
    """Casimir force in the high temperature limit given by the PFA.

    Parameters:
        LbyR: L/R

    Return value:
        F_PFA/(kb*T/R)
    """
    return -zeta3/(8*LbyR**2)


def E_HT_PFA(LbyR):
    """Casimir interaction energy in the high temperature limit given by the
    PFA.

    Parameters:
        LbyR: L/R

    Return value:
        F_PFA/(kb*T)
    """
    return -zeta3/(8*LbyR)


def F_HT_cut(eps, a=1e-4, b=1):
    """Compute LbyR so that
         1-F_exact/F_PFA = eps.
    The function uses the bisection method from scipy.

    Parameters:
        eps: accuracy
        a: left border for bisection
        b: right border for bisection

    Return value:
        L/R
    """
    ratio = lambda x: 1-F_HT_exact(x)/F_HT_PFA(x)-eps
    return bisect(ratio, a, b)


def E_HT_cut(eps, a=1e-4, b=1):
    """Compute LbyR so that
         1-E_exact/E_PFA = eps.
    The function uses the bisection method from scipy.

    Parameters:
        eps: accuracy
        a: left border for bisection
        b: right border for bisection

    Return value:
        L/R
    """
    ratio = lambda x: 1-E_HT_exact(x)/E_HT_PFA(x)-eps
    return bisect(ratio, a, b)

if __name__ == "__main__":
    print("force 1%:  ", 1/F_HT_cut(0.01))
    print("force 0.5%:", 1/F_HT_cut(0.005))
    print("force 0.3%:", 1/F_HT_cut(0.003))

    print()

    print("energy 2%:  ", 1/E_HT_cut(0.02))
    print("energy 1%:  ", 1/E_HT_cut(0.01))
    print("energy 0.5%:", 1/E_HT_cut(0.005))
