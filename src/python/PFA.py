import numpy as np
from math import exp, fsum, pi, sqrt, log1p, isinf
from scipy.integrate import quad

# constants
c       = 299792458       # speed of light [m/s]
kB      = 1.38064852e-23  # Boltzmann constant [m² kg / (s² K)]
hbar    = 1.0545718e-34   # hbar [J s]
hbar_eV = 6.582119514e-16 # hbar [eV s / rad]

def psum(f, L, T, epsrel=1e-10):
    """Evaluate the primed sum
         k_B T \sum_n^\prime f(L,xi_n)
    where
        xi_n = 2*\pi*k_B*T/\hbar
    are the Matsubara frequencies and f is a function that depends on the
    separation L and the Matsubara frequency xi_n.
    
    If T=0, the sum becomes an integral. This function also evaluates T=0
    correctly.
    """
    if T == 0:
        # calculate integral from 0 to inf of hbar*f(t)/(2*pi). 
        # scale integrand
        alpha = c/L
        integrand = lambda x: f(L, alpha*x)
        I,err = quad(integrand, 0, np.inf, epsabs=0, epsrel=epsrel)
        return alpha*hbar*I/(2*pi)
    else:
        terms = []
        n = 0
        while True:
            xi_n = 2*pi*n*kB*T/hbar
            I = f(L, xi_n)
            terms.append(I)
            n += 1

            if abs(I/terms[0]) < epsrel:
                terms[0] /= 2
                return kB*T*fsum(terms)



class PFA:
    def __init__(self, R, T, epsm1):
        """Initialize PFA object
        R:     radius of sphere in m
        T:     temperature in Kelvin
        epsm1: dielectric function minus 1, epsm1(xi) = eps-1
        """
        self.R = R
        self.T = T
        self.epsm1 = epsm1


    def rp(self, xi, kappa):
        """Fresnel coefficients
        Compute Fresnel coefficients r_TE and r_TM.
        
        If xi=0, r_TE=0 and r_TM=1 will be returned (Drude).
        If epsm1(xi) is inf, r_TE=-1 and r_TM=1 will be returned for xi > 0.

        Returns r_TE,r_TM
        """
        if xi == 0:
            return 0,1 # Drude

        epsm1 = self.epsm1(xi)
        if isinf(epsm1): # PR
            return -1, 1

        beta = sqrt(1+ (xi/(c*kappa))**2*epsm1 )

        rTE = (1-beta)/(1+beta)
        rTM = (epsm1+1-beta)/(epsm1+1+beta)

        return rTE, rTM


    def dF_xi(self, L, xi, epsabs=0, epsrel=1e-10):
        def f(t,xi):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            rTE2, rTM2 = rTE**2, rTM**2
            exp_f = exp(-2*t)
            return t**2*exp_f*( rTE2/(1-rTE2*exp_f) + rTM2/(1-rTM2*exp_f) )

        I,err = quad(f, xi*L/c, np.inf, args=(xi,), epsabs=epsabs, epsrel=epsrel)
        return I


    def dF(self, L):
        """ dF = dF/dL = - d²F/dL²"""
        return self.R/L**3*2*psum(self.dF_xi, L, self.T)


    def F_xi(self, L, xi, epsabs=0, epsrel=1e-10):
        def f(t,xi):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            exp_f = exp(-2*t)
            return t*(log1p(-rTE**2*exp_f)+log1p(-rTM**2*exp_f))

        I,err = quad(f, xi*L/c, np.inf, args=(xi,), epsabs=epsabs, epsrel=epsrel)
        return I


    def F(self, L):
        """Casimir force, F = -dE/dL"""
        return self.R/L**2*psum(self.F_xi, L, self.T)


    def E_xi(self, L, xi, epsabs=0, epsrel=1e-10):
        def f(t,xi):
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            exp_f = exp(-2*t)
            return (log1p(-rTE**2*exp_f)+log1p(-rTM**2*exp_f))*(t-xi*L/c)

        I,err = quad(f, xi*L/c, np.inf, args=(xi,), epsabs=epsabs, epsrel=epsrel)
        return I


    def E(self,L):
        """Casimir free energy E"""
        return self.R/L*psum(self.E_xi, L, self.T)
