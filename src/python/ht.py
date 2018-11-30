import itertools
from functools import lru_cache
from math import log,exp,lgamma,sqrt,fsum
import numpy as np
from scipy.integrate import quad
from scipy.special import ive

class Plasma:
    # Implementation of the high-temperature case for plasma. this serves as a
    # reference implementation.

    def __init__(self, L, R, omegap=9):
        """
        Compute the Casimir free energy in the high-temperature limit for the
        plane-sphere geometry. Supports Drude and plasma model.

        parameters:
            L: separation between plane and sphere
            R: radius of sphere
            omegap: plasma frequency in eV
        """
        hbar_eV = 6.582119514e-16 # hbar [eV s / rad]
        c       = 299792458       # speed of light [m/s]

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
            f = sqrt(x**2+self.beta**2)
            return exp(nu*log(x)-x-lgamma(1+nu)) * (f-x)/(f+x)
            
        I1, err1 = quad(integrand, 0, nu)
        I2, err2 = quad(integrand, nu, float("inf"))
        return I1+I2


    def elem(self, l1,l2,m):
        """Matrix elements of the round-trip operator for the plasma model"""
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
        """logdet(Id-M) for plasma model"""
        D = np.eye(ldim) - self.M(m,ldim)
        s,v = np.linalg.slogdet(D)
        assert(s == 1)
        return v


    def E(self, verbose=False, cutoff=1e-8, eta=8):
        """MM contribution to the free energy in the high-temperature limit for
        the plasma model in units of kb*T"""
        ldim = int(eta/self.LbyR)

        terms = []
        for m in itertools.count():
            v = self.logdetD(m,ldim)
            terms.append(v)

            if verbose:
                print("#", m, v)

            if v/terms[0] < cutoff:
                terms[0] /= 2
                return fsum(terms)
 

def __dirichlet(L,R,lmax=10000,lmin=0):
    x = L/R
    Z = 1+x-np.sqrt(x*(x+2))
    l = np.arange(lmin,lmax+1)

    # dZ/dL
    dZ = (1-(x+1)/np.sqrt(x*(x+2)))/R

    # d²Z/dL²
    dZ2 = np.sqrt(x*(x+2))/(4*x**2+4*x**3+x**4)/R**2

    A = Z**(2*l) # Z^(2l)
    q  = 2*l+1   # 2l+1
    q2 = q**2    # (2l+1)²

    # E = kb T/2 Σ (2l+1) log(1-Z^(2l+1))
    E = np.sum(q*np.log1p(-A*Z))/2

    # B = Z^(2l)/(1-Z^(2l+1))
    B = A/(1-A*Z)

    # F = -dE/dL
    F = np.sum(q2*B)*dZ/2

    # dF = F' = dF/dL = -d²E/dL²
    dF = np.sum(q2*B*( dZ**2*(2*l/Z + q*B) + dZ2 ))/2

    return E,F,dF


def dirichlet(L,R,lmax=10000):
    """Computes the Casimir interaction (i.e., the free energy E, the force
    F=-dE/dL, and the force gradient F'=-d²E/dL²) for the plane-sphere geometry
    in the high-temperature limit for Dirichlet boundary conditions.

    The function evaluates Eq. (3) of Bimonte, Emig, PRL 109, 160403 (2012).

    Parameters:
        L: separation between plane and sphere
        R: radius of sphere
        lmax: truncation of sum

    Returns:
        (E, F, F')
    """
    return __dirichlet(L,R,lmax=lmax,lmin=0)


def drude(L,R,lmax=10000):
    """Computes the Casimir interaction (i.e., the free energy E, the force
    F=-dE/dL, and the force gradient F'=-d²E/dL²) for the plane-sphere geometry
    in the high-temperature limit for Drude boundary conditions.

    The function evaluates Eq. (8) of Bimonte, Emig, PRL 109, 160403 (2012).

    Parameters:
        L: separation between plane and sphere
        R: radius of sphere
        lmax: truncation of sum

    Returns:
        (E, F, F')
    """
    # E1 = kb T/2 Σ (2l+1) log(1-Z^(2l+1))  where l=1,2,..,lmax
    E1, F1, dF1 = __dirichlet(L,R,lmax=lmax,lmin=1)

    x = L/R
    Z = 1+x-np.sqrt(x*(x+2))
    l = np.arange(1,lmax+1)

    # dZ/dL
    dZ = (1-(x+1)/np.sqrt(x*(x+2)))/R

    # d²Z/dL²
    dZ2 = np.sqrt(x*(x+2))/(4*x**2+4*x**3+x**4)/R**2

    A = Z**(2*l) # Z^(2l)
    q  = 2*l+1   # 2l+1
    q2 = q**2    # (2l+1)²

    # B = Z^(2l)/(1-Z^(2l+1))
    B = A/(1-A*Z)

    # E2 = kb T/2 log(1-(1-Z²)Σ Z^(2l+1) (1-Z^(2l))/(1-Z^(2l+1)))
    S = np.sum((1-A)*Z*B)
    E2 = np.log1p(-(1-Z**2)*S)/2

    # dS/dZ
    dS = np.sum(B*((1-A)*q+Z*(1-A)*q*B-2*A*l))
    
    # d²S/dZ²
    dS2 = np.sum( 2*B*(B*(1-A)*q*(3*l+1+q*Z*B) - 4*A*l**2/Z - A/Z*l*q + 1/Z*(1-A)*l*q - 2*B*A*l*q))

    # F2 = -dE2/dL
    F2 = -dZ/2*(2*Z*S-(1-Z**2)*dS)/(1-(1-Z**2)*S)

    # dF2 = -d²E2/dL²
    dF2 = -(dZ2*(2*S*Z-(1-Z**2)*dS)/(1-(1-Z**2)*S) + dZ**2*( (1-(1-Z**2)*S)*(4*dS*Z+2*S-(1-Z**2)*dS2) - (2*S*Z-(1-Z**2)*dS)*(2*S*Z-(1-Z**2)*dS) )/(1-(1-Z**2)*S)**2)/2

    return E1+E2, F1+F2, dF1+dF2


def logdet_MM_0(x):
    """Compute logdet(Id-M^0_MM) for xi=0 and m=0"""
    l = np.arange(20000)
    Z = 1/(1+x+np.sqrt(x*(2+x)))
    return np.sum(np.log1p(-Z**(2*l+3)))


def logdet_EE_0(x):
    """Compute logdet(Id-M^0_EE) for xi=0 and m=0"""
    l = np.arange(1,20000)
    Z = 1/(1+x+np.sqrt(x*(2+x)))
    sum1 = np.sum(np.log1p(-Z**(2*l+1)))
    sum2 = np.sum(Z**(2*l+1)*(1-Z**(2*l))/(1-Z**(2*l+1)))

    return sum1 + np.log1p(-(1-Z**2)*sum2)


if __name__ == "__main__":
    L = 1
    R = 20

    ht = Plasma(L,R)
    plasma = ht.E(verbose=True, cutoff=1e-12)

    print(plasma)
