import ctypes
from math import pi
from scipy.integrate import dblquad
from sys import path
path.append("..")

from PFA import PFA, hbarc, c, hbar_eV

def theta_teo(L, R, omegap=9, gamma=0.035, epsrel=1e-11):
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
        epsrel: relative accuracy for integration

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

    I, err = dblquad(integrand, 0, float("inf"), lambda x: x, lambda x: float("inf"), args=(L,R,omegap,gamma), epsrel=epsrel)
    f = R/L*(1+R/L)/(4*pi)
    E1 = -I*f

    theta = R/L*E1/E0

    return theta


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute linear correction to PFA in the plane-sphere geometry for T=0")
    parser.add_argument("--omegap", type=float, default=9, help="plasma frequency in eV")
    parser.add_argument("--gamma", type=float, default=0, help="relaxation frequency in eV")
    parser.add_argument("-L", type=float, required=True, help="separation between sphere and plate in m")
    parser.add_argument("-R", type=float, required=True, help="radius of sphere in m")
    args = parser.parse_args()

    print("L      = %g" % args.L)
    print("R      = %g" % args.R)
    print("omegap = %g" % args.omegap)
    print("gamma  = %g" % args.gamma)

    theta = theta_teo(args.L, args.R, args.omegap, args.gamma)
    print("theta  = %.15g" % theta)
