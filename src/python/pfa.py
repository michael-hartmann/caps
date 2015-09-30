#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
from math import pi,exp,isinf
from scipy.integrate import quad
from scipy.special import zetac

def integrand(x, LbyR, T, prec=1e-15):
    """determine integrand for PFA expression (7.6) in [1]

    Parameters:
    -----------
    x : float, integration variable
    LbyR : float, ratio of surface-surface distance to sphere radius
    T : float, temperature scaled according to (5.41) in [1]
    prec : float, precision used in convergence criterion of the sum

    Returns:
    --------
    integrand including the factor 1/(L/R) in the prefactor

    Notes:
    ------
    In order to avoid the evaluation of the di- and trilogarithm, the
    integrand has been resummed. In the Matsubara sum, the n=0 term is
    accounted for explicitly by means of Li_3(1)/2 = ζ(3)/2. Expressing
    the tri- and dilogarithm in terms of their series expansion (cf.
    (A.23) in [1]), the Matsubara sum yields a geometric sum and its
    derivative, respectively, for which a closed expression can be given.
    It then remains one sum to be done.

    References:
    -----------
    [1] M. Hartmann, master thesis (Univ. Augsburg, 2014)
        https://github.com/michael-hartmann/libcasimir/raw/master/thesis.pdf

    """
    n = 1
    sum = 0.5*(1+zetac(3))
    alpha1 = 2*T*x*LbyR/(1+LbyR)
    while True:
        alpha = alpha1*n
        arg = exp(-alpha)
        value = arg/(1-arg)*(1/n**3+alpha1/(n**2*(1-arg)))
        sum += value
        if value/sum < prec:
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
    from sys import argv, stderr, stdout, exit

    def usage(self, stream=stdout):
        print("%s LbyR  T" % (self),         file=stream)
        print("\tLbyR: ratio L/R, LbyR > 0", file=stream)
        print("\tT:    temperature, T > 0",  file=stream)

    if len(argv) < 3:
        usage(argv[0])
        exit(1)

    try:
        LbyR = float(argv[1])
        T    = float(argv[2])

        if LbyR < 0:
            raise BaseException("LbyR must be positive")
        if T < 0:
            raise BaseException("T must be positive")
    except:
        usage(argv[0], stream=stderr)
        exit(1)

    print("%+.15g" % pfa(LbyR, T))
