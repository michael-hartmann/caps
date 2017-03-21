import numpy as np
from scipy import interpolate

hbar_eV = 6.582119514e-16 # eV/s

def get_epsilonm1(filename, omega_p=9.0, gamma=0.03):
    data = np.loadtxt(filename)

    data_xi  = data[:,0]
    data_eps = data[:,1]
    xi_min, xi_max = data_xi[0], data_xi[-1]

    f = interpolate.interp1d(data_xi, data_eps, kind="linear")

    omega_p /= hbar_eV
    gamma   /= hbar_eV

    def epsilonm1(xi):
        if xi < xi_min:
            return omega_p**2/(xi*(xi+gamma))
        
        if xi > xi_max:
            return 0

        return f(xi)-1

    return epsilonm1
