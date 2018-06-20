import ctypes
from math import sqrt, pi
from itertools import count
from scipy.integrate import dblquad

clib = ctypes.CDLL("./integrand.so")
clib.integrand.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
clib.integrand.restype = ctypes.c_double

c       = 299792458       # speed of light [m/s]
hbar_eV = 6.582119514e-16 # hbar [eV s / rad]

def epsilon(xi, omegap, gamma):
    xi_eV = xi*hbar_eV
    return 1+omegap**2/(xi_eV*(xi_eV+gamma))

def integrand(x, xi_, L, R, omegap, gamma):
    assert x >= xi_
    e = L/R # L/R

    xi = xi_*c/L
    eps = epsilon(xi, omegap, gamma)

    return clib.integrand(xi_, x, eps, e)

omegap = 9
gamma = 0 #0.035
R = 150e-6
L = 300e-9

#def f(x,y,L,R,omegap,gamma):
#    l = y
#    tau = 2*e*l/x
#    dx = 2*e/tau
#    return dx*integrand(tau, l, L, R, omegap, gamma)

I, err = dblquad(integrand, 0, float("inf"), lambda x: x, lambda x: float("inf"), args=(L,R,omegap,gamma), epsrel=1e-8)

f = R/L*(1+R/L)/(4*pi)
E0 = -9036.32256859 #-31147.0818086
E1 = -I*f
#print("pfa", pfa)

print(omegap/hbar_eV*L/c)
theta = R/L*E1/E0
print(theta)
