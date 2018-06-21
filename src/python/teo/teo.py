import ctypes
from math import pi
from scipy.integrate import dblquad
from sys import path
path.append("..")
from uncertainties import ufloat

from PFA import PFA, hbarc, c, hbar_eV

# This script computes the linear correction Θ_1 at T=0. The correction looks
# like
#       E = E_PFA (1+L/R*Θ_1)
# where Θ_1 is a function of L but independent of R.
#
# To compute the linear correction, we use the formulas implemented in Ref [1].
# The Casimir energy is given as
#       E = E0 + E1
# where E0 corresponds to the PFA result, E0=E_PFA, and E1 is the correction.
# The linear term can be computed by means of
#       Θ_1 = R/L * E1/E0.
# The linear correction Θ_1 is independent and one gets the same value for Θ_1
# independent of the choice of R. However, the integration error gets
# multiplied by R/L. For this reason high aspect ratios increase the
# uncertainty of Θ_1. On the other hand, the integration has problems when the
# aspect ratio is not sufficiently high.
#
# The actual integrand contains a lot of computations and an infinite sum. For
# performance reasons, it is implemented in a C library, see integrand.c. The
# library is called using ctypes.
#
# [1] Teo, Material dependence of Casimir interaction between a sphere and a
# plate: First analytic correction beyond proximity force approximation,
# https://doi.org/10.1103/PhysRevD.88.045019

def theta_teo(L, R, omegap=9, gamma=0.035, epsrel=1e-11, epsabs=1e-8):
    """Compute linear correction Θ_1 for T=0

    The linear correction is given as
        E ≈ E_0 + E_1 = E_0*(1+E_1/E_0) = E_0*(1+Θ_1*x)
    where x=L/R and E_0 corresponds to the PFA result.

    Parameters
    ----------
        L:      separation between sphere and plate in m
        R:      radius of the sphere in m
        omegap: plasma frequency in eV
        gamma:  relaxation frequency in eV
        epsrel: relative tolerance for integration
        epsabs: absolute tolerance for integration

    Returns
    -------
    Θ_1: linear correction
    """

    # open library
    clib = ctypes.CDLL("./integrand.so")
    clib.integrand.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    clib.integrand.restype = ctypes.c_double

    epsm1 = lambda xi: (omegap/hbar_eV)**2/(xi*(xi+gamma/hbar_eV))
    pfa = PFA(R, 0, epsm1)
    E0 = pfa.E(L)*(L+R)/hbarc

    omegap_rads = omegap/hbar_eV
    gamma_rads  = gamma/hbar_eV

    def integrand(x, xi_, L, R, omegap, gamma):
        xi = xi_*c/L
        eps = 1+omegap_rads**2/(xi*(xi+gamma_rads))
        return clib.integrand(xi_, x, eps, L/R)

    I, err = dblquad(integrand, 0, float("inf"), lambda x: x, lambda x: float("inf"), args=(L,R,omegap,gamma), epsrel=epsrel, epsabs=epsabs)
    I = ufloat(I,err)
    f = R/L*(1+R/L)/(4*pi)
    E1 = -I*f

    theta = R/L*E1/E0

    return theta


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute linear correction to PFA in the plane-sphere geometry for T=0")
    parser.add_argument("--omegap", type=float, default=9, help="plasma frequency in eV")
    parser.add_argument("--gamma", type=float, default=0, help="relaxation frequency in eV")
    parser.add_argument("--epsrel", type=float, default=1e-11, help="relative tolerance for integration")
    parser.add_argument("--epsabs", type=float, default=1e-8, help="absolute tolerance for integration")
    parser.add_argument("-L", type=float, required=True, help="separation between sphere and plate in m")
    parser.add_argument("-R", type=float, required=True, help="radius of sphere in m")
    args = parser.parse_args()

    print("L      = %g" % args.L)
    print("R      = %g" % args.R)
    print("omegap = %g" % args.omegap)
    print("gamma  = %g" % args.gamma)

    theta = theta_teo(args.L, args.R, omegap=args.omegap, gamma=args.gamma, epsrel=args.epsrel, epsabs=args.epsabs)
    print("theta  = ", theta)
