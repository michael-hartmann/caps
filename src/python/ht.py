import itertools
from functools import lru_cache
from math import *
import numpy as np
from scipy.integrate import quad
from scipy.special import ive

hbar_eV = 6.582119514e-16 # hbar [eV s / rad]
c       = 299792458       # speed of light [m/s]

class CasimirHT:
    def __init__(self, L, R, omegap=9):
        """
        Compute the Casimir free energy in the high-temperature limit for the
        plane-sphere geometry. Supports Drude and plasma model.

        parameters:
        L: separation between plane and sphere
        R: radius of sphere
        omegap: plasma frequency in eV
        """
        self.LbyR  = L/R                        # L/R
        self.y     = log(R/(L+R)/2)             # log(R/(2(L+R))
        self.alpha = omegap/(hbar_eV*c)*R       # α=ωp*R/c
        self.beta  = 2*omegap/(hbar_eV*c)*(R+L) # β=2*ωp*(L+R)/c

    @lru_cache(maxsize=32000)
    def ratio(self, l):
        """Compute the ratio I_{l+1/2}(alpha) / I_{l-1/2}(alpha)"""
        return ive(l+0.5,self.alpha)/ive(l-0.5,self.alpha)


    @lru_cache(maxsize=64000)
    def integral(self, nu):
        """
        Compute the integral 
            I = x^nu*exp(-x)*sqrt(1+(β/2)²-1)/sqrt(1+(β/2)²+1)/nu!
        """
        def integrand(x):
            f = sqrt(1+(self.beta/x)**2)
            return exp(nu*log(x)-x-lgamma(1+nu)) * (f-1)/(f+1)
            
        I1, err1 = quad(integrand, 0, nu)
        I2, err2 = quad(integrand, nu, float("inf"))
        return I1+I2


    def elem(self, l1,l2,m):
        """Matrix elements of the round-trip operator for plasma"""
        f1 = sqrt(l1/(l1+1)*l2/(l2+1))
        f2 = sqrt(1-(2*l1+1)/self.alpha*self.ratio(l1))
        f3 = sqrt(1-(2*l2+1)/self.alpha*self.ratio(l2))

        return f1*f2*f3* self.integral(l1+l2) * exp((l1+l2+1)*self.y + lgamma(1+l1+l2)-(lgamma(1+l1+m)+lgamma(1+l1-m)+lgamma(1+l2+m)+lgamma(1+l2-m))/2)


    def M(self,m,ldim):
        """Compute round-trip matrix for plasma model"""
        A = np.zeros((ldim,ldim))
        offset = max(1,m)
        for l1 in range(offset,ldim+1):
            for l2 in range(l1,ldim+1):
                A[l1-offset,l2-offset] = A[l2-offset,l1-offset] = self.elem(l1,l2,m)

        return A


    def logdetD(self,m,ldim):
        D = np.eye(ldim) - self.M(m,ldim)
        s,v = np.linalg.slogdet(D)
        assert(s == 1)
        return v


    def plasma(self, verbose=False, cutoff=1e-8):
        terms = []
        for m in itertools.count():
            v = self.logdetD(m,ldim)
            terms.append(v)

            if verbose:
                print("#", m, v)

            if v/terms[0] < cutoff:
                terms[0] /= 2
                return fsum(terms)
 

    def drude(self, lmax=10000):
        LbyR = self.LbyR
        Z = (1+LbyR)*(1-sqrt(LbyR*(LbyR+2)/(LbyR**2+2*LbyR+1)))

        l = np.arange(1,lmax+1)
        term1 = np.sum((2*l+1)*np.log1p(-Z**(2*l+1)))
        term2 = log1p(-(1-Z**2)*np.sum(Z**(2*l+1)*(1-Z**(2*l))/(1-Z**(2*l+1))))

        return (term1+term2)/2


if __name__ == "__main__":
    omegap = 9
    R = 10e-6
    L = 200e-9
    LbyR = L/R
    ldim = int(1+8/LbyR)

    ht = CasimirHT(L,R,omegap=omegap)
    plasma = ht.plasma(verbose=True)
    drude  = ht.drude()

    print()
    print(drude, plasma, drude+plasma)
