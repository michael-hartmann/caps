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

def integrand(tau, l, L, R, omegap, gamma):
    e = L/R # L/R
    xi = c/R*l*sqrt(1-tau**2)/tau # xi in rad/s
    eps = epsilon(xi, omegap, gamma)

    return clib.integrand(tau, l, eps, e)

omegap = 9
gamma = 0.035
R = 150e-6
L = 150e-9

#def f(x,y,L,R,omegap,gamma):
#    l = y
#    tau = 2*e*l/x
#    dx = 2*e/tau
#    return dx*integrand(tau, l, L, R, omegap, gamma)

I, err = dblquad(integrand, 0, float("inf"), lambda x: 0, lambda x: 1, args=(L,R,omegap,gamma), epsrel=1e-6)
f = (1+L/R)/(4*pi)
pfa = -I*f
err *= f
print("I", I)
print("pfa", pfa)
print("err", abs(err/pfa))
