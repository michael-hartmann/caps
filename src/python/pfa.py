#!/usr/bin/python

from __future__ import division
from math import *
from scipy.integrate import quad
from scipy.special import zetac

def integrand(x, LbyR, T):
    """alternative integrand for PFA and finite temperature

    """
    n = 1
    sum = 0.5*(1+zetac(3))
    alpha1 = 2*T*x*LbyR/(1+LbyR)
    while True:
        alpha = alpha1*n
        arg = exp(-alpha)
        value = arg/(1-arg)*(1/n**3+alpha1/(n**2*(1-arg)))
        sum += value
        if value/sum < 1e-15:
            return sum/(x**2*LbyR)
        n += 1

def pfa(LbyR, T):
    """Calculate free energy according to PFA for perfect reflectors for L/R
    and T"""
    if isinf(T) and T > 0:
        return -1.2020569031595942/4/(LbyR+LbyR**2)
    if T > 0:
        I = quad(integrand, 1, 1+1/LbyR, args=(LbyR, T), epsrel=1e-12)
        return -T/(4*pi)*I[0]
    elif T == 0:
        # T = 0
        return -pi**3/720*(1+2*LbyR)/(LbyR**2+LbyR**3)
    else:
        raise BaseException("invalid value for T")


if __name__ == "__main__":
    from sys import argv, exit

    def usage(self):
        print "%s LbyR, T" % (self)
        print "\tLbyR: ratio L/R, LbyR > 0"
        print "\tT:    temperature, T > 0"

    if argv < 3:
        usage(argv[0])
        exit(1)

    try:
        LbyR = float(argv[1])
        T    = float(argv[2])

        if LbyR < 0:
            raise BaseException("LbyR < 0")
        if T < 0:
            raise BaseException("R < 0")
    except:
        usage(argv[0])
        exit(1)

    print pfa(LbyR, T)
