import ctypes
from math import sqrt, pi
from scipy.integrate import dblquad
from sys import path
path.append("..")

from PFA import PFA, hbarc, c, hbar_eV

clib = ctypes.CDLL("./integrand.so")
clib.integrand.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
clib.integrand.restype = ctypes.c_double

def integrand(x, xi_, L, R, omegap, gamma):
    xi = xi_*c/L

    xi_eV = xi*hbar_eV
    eps = 1+omegap**2/(xi_eV*(xi_eV+gamma))

    return clib.integrand(xi_, x, eps, L/R)

def theta_teo(L,R,omegap=9,gamma=0.035):
    epsm1 = lambda xi: (omegap/hbar_eV)**2/(xi*(xi+gamma/hbar_eV))
    pfa = PFA(R, 0, epsm1)
    E0 = pfa.E(L)*(L+R)/hbarc

    I, err = dblquad(integrand, 0, float("inf"), lambda x: x, lambda x: float("inf"), args=(L,R,omegap,gamma), epsrel=1e-11)
    f = R/L*(1+R/L)/(4*pi)
    E1 = -I*f

    theta = R/L*E1/E0

    return theta


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute PFA in the plane-sphere geometry")
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
