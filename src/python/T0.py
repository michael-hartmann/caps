import numpy as np
from scipy.integrate import quad
from argparse import ArgumentParser
from libcasimir import Casimir
from math import fsum,ceil
from sys import stderr

def integrand(xi, LbyR, ldim, threshold=None, precision=1e-10):
    terms = []
    m = 0
    while True:
        casimir = Casimir(LbyR, ldim=ldim)
        if threshold != None:
            casimir.set_threshold(threshold)
        logdetD_m = casimir.logdetD(xi, m)
        del casimir

        terms.append(logdetD_m)
        if m == 0 and terms[m] == 0:
            break
        elif abs(terms[m]/terms[0]) < precision:
            break
        m = m+1

    terms[0] = terms[0]/2
    terms.sort()
    value = fsum(terms)

    print("# xi*(L+R)/c=%.15g, logdetD(xi)=%.15g" % (xi, value))
    return value


if __name__ == "__main__":
    parser = ArgumentParser(description="Calculate free Casimir energy for T=0")
    parser.add_argument("-x", "--LbyR", action="store", dest="LbyR", type=float,   help="aspect ratio x=L/R")
    parser.add_argument("-L", "--ldim", action="store", dest="ldim", type=int)

    parser.add_argument("--eta",    action="store", dest="eta",  type=float, default=6)
    parser.add_argument("--lmin",   action="store", dest="lmin", type=float, default=30)
    parser.add_argument("--epsrel", action="store", dest="epsrel", type=float, default=1e-6)

    parser.add_argument("-p", "--precision",  action="store", dest="precision", type=float, default=1e-10)
    parser.add_argument("-t", "--threshold" , action="store", dest="threshold", type=float)

    args = parser.parse_args()

    LbyR = args.LbyR
    if LbyR <= 0:
        print("-x, --LbyR: argument must be positive")
        exit(1)

    if args.ldim != None:
        if args.ldim < 1:
            print("-L, --ldim: argument must be positive integer", file=stderr)
            exit(1)
        ldim = args.ldim
    else:
        if args.lmin < 1:
            print("--lmin: argument must be a positive integer", file=stderr)
            exit(1)
        if args.eta <= 0:
            print("--eta: argument must be positive", file=stderr)
            exit(1)
        ldim = int(ceil(max(args.lmin, args.eta/LbyR)))

    if args.precision != None:
        if args.precision <= 0:
            print("-p, --precision: argument must be positive")
            exit(1)
    precision = args.precision

    if args.epsrel < 0:
        print("--epsrel: argument must be positive")
        exit(1)
    epsrel = args.epsrel

    if args.threshold != None:
        if threshold <= 0:
            print("-t, --threshold: argument must be positive")
            exit(1)
    threshold = args.threshold

    print("#", args)
    print("#")

    integral,abserr,infodict = quad(integrand, 0, np.inf, epsabs=0, epsrel=epsrel, full_output=True, args=(LbyR,ldim,threshold,precision))
    F = integral/np.pi

    print("#")
    print("# neval = %d" % infodict["neval"])
    print("#")
    print("# L/R, ldim, F(T=0)*(L+R)/(hbar*c), error (absolute), error (relative)")
    print("%.15g, %d, %.15g, %g, %g" % (LbyR, ldim, F, abserr/np.pi, abserr/integral))
