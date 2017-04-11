import functools
import numpy as np
from scipy.integrate import quad
from scipy.special import ive
from math import *

hbar    = 1.0545718e-34   # m² kg / s
hbar_eV = 6.582119514e-16 # eV s
kB      = 1.38064852e-23  # m² kg / ( K s² )
c       = 299792458       # m/s


R = 150e-6
L =  30e-6
omega_p = 9/hbar_eV
lambda_p = 2*pi*c/omega_p

alpha = 2*pi*R/lambda_p

x = R/(R+L)
Ltilde = 1/x

@functools.lru_cache(maxsize=100, typed=False)
def integral(nu):
    def integrand(t):
        beta = sqrt(1+(2*alpha*Ltilde/t)**2)
        return t**nu*(beta-1)/(beta+1)*exp(-t)

    I,err = quad(integrand, 0, np.inf, epsrel=1e-12)
    return I


m = 1
lmin = 1
lmax = 30

dim = lmax-lmin+1
M = np.zeros((dim,dim))

for l1 in range(lmin,lmax+1):
    ratio1  = ive(l1+0.5,alpha)/ive(l1-0.5,alpha)
    factor1 = 1-(2*l1+1)/alpha*ratio1

    for l2 in range(lmin,lmax+1):

        ratio2  = ive(l2+0.5,alpha)/ive(l2-0.5,alpha)
        factor2 = 1-(2*l2+1)/alpha*ratio2

        factor = sqrt(factor1*factor2)

        nu = l1+l2
        
        #elem = 1/(2*Ltilde)**(l1+l2+1)*l1/(l1+1)/sqrt(factorial(l1-m)*factorial(l1+m)*factorial(l2-m)*factorial(l2+m)) * factor*integral(nu)
        elem = 1/(2*Ltilde)**(l1+l2+1)*sqrt(l1*l2/((l1+1)*(l2+1)))/sqrt(factorial(l1-m)*factorial(l1+m)*factorial(l2-m)*factorial(l2+m)) * factor*integral(nu)

        M[l1-lmin,l2-lmin] = elem

D = np.eye(dim)-M

print(np.linalg.slogdet(D))
