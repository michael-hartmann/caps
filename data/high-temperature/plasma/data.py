import numpy as np
from sys import argv,path
from glob import glob
from math import sqrt, exp 
import scipy.special
from scipy.integrate import quad

c       = 299792458       # speed of light [m/s]
hbar_eV = 6.582119514e-16 # hbar [eV s / rad]

# dilogarithm
Li2 = lambda x: scipy.special.spence(1-x)

def E_plasma_ht(L,R,omegap_eV):
    omegap = omegap_eV/hbar_eV
    def f(t):
        kappa = t/L
        beta = sqrt(1 + (omegap/(c*kappa))**2)
        rTM, rTE = 1, (1-beta)/(1+beta)
        exp_f = exp(-2*t)
        return Li2(rTE**2*exp_f)+Li2(rTM**2*exp_f)

    I,err = quad(f, 0, float("inf"), epsabs=0, epsrel=1e-12)
    return -0.5*R/L*I/2


data = []
for filename in glob(argv[1]+"/raw-*.out"):
    LbyR, L, R, ldim, omegap, F_drude, F_PR, F_plasma = np.loadtxt(filename, delimiter=",")
    data.append((1/LbyR, F_drude, F_PR, F_plasma))

print("# high-temperature limit")
print("# R: %g m" % R)
print("# omegap: %g eV" % omegap)
print("#")
print("# R/L, F_Drude/(kB*T), F_PR/(kB*T), F_plasma/(kB*T), F_plasma/F_PFA_plasma")

for RbyL, drude, pr, plasma in sorted(data):
    pfa = E_plasma_ht(R/RbyL,R,omegap)
    print("%.8g, %.14g, %.14g, %.14g, %.14g" % (RbyL, drude, pr, plasma, plasma/pfa))
