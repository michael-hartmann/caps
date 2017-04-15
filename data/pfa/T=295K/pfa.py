import numpy as np
from math import exp, fsum, pi, sqrt
from scipy.integrate import quad
from scipy import interpolate

# constants
c    = 299792458       # speed of light [m/s]
kB   = 1.38064852e-23  # Boltzmann constant [m² kg / (s² K)]
hbar = 1.0545718e-34   # hbar [J s]
hbar_eV = 6.582119514e-16 # hbar [eV s]

class PFA:
    def __init__(self, R, T, epsilonm1, omegap=0, args=None):
        """R in m, T in K, omegap and gamma in eV"""
        self.omegap = omegap
        self.R = R
        self.T = T
        self.epsilonm1 = epsilonm1
        self.args = args


    def rp(self, xi, kappa):
        if xi == 0:
            # plasma    
            if self.omegap:
                beta = sqrt(1+(self.omegap/(c*kappa))**2)
                rTE = (1-beta)/(1+beta)
                return rTE,1
            else:
                # Drude
                return 0,1

        epsm1 = self.epsilonm1(xi, self.args)
        beta = sqrt(1+ (xi/(c*kappa))**2*epsm1 )

        rTE = (1-beta)/(1+beta)
        rTM = (epsm1+1-beta)/(epsm1+1+beta)

        return rTE, rTM


    def P(self, L):
        def f(t,xi):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            rTE2, rTM2 = rTE**2, rTM**2
            exp_f = exp(-2*t)
            return t**2*exp_f*( rTE2/(1-rTE2*exp_f) + rTM2/(1-rTM2*exp_f) )

        n = 0
        terms = []
        while True:
            xi_n = 2*pi*n*kB*self.T/hbar
            I,err = quad(f, xi_n*L/c, np.inf, args=(xi_n,), epsabs=0, epsrel=1e-12)
            terms.append(I)

            if abs(I/terms[0]) < 1e-10:
                break

            n += 1

        terms[0] /= 2
        return -(kB*self.T)/(2*pi)/L**3*2*fsum(terms)


def get_epsilonm1(filename, omega_p, gamma):
    data = np.loadtxt(filename)

    data_xi  = data[:,0]
    data_eps = data[:,1]
    xi_min, xi_max = data_xi[0], data_xi[-1]

    f = interpolate.interp1d(data_xi, data_eps, kind="linear")

    def epsilonm1(xi,args):
        if xi < xi_min:
            return omega_p**2/(xi*(xi+gamma))
     
        if xi > xi_max:
            # XXX
            return 0

        return f(xi)-1

    return epsilonm1



def pressure(L,R,T):
    omegap = 9/hbar_eV    # plasma frequency 9eV
    gamma  = 0.03/hbar_eV # dissipation 0.03
    filename = "../../materials/GoldEpsIm.dat"

    pfa = PFA(R,T,get_epsilonm1(filename, omegap, gamma))
    P_drude = pfa.P(L)*1000 # in mPa
    
    pfa = PFA(R,T,get_epsilonm1(filename, omegap, gamma), omegap=omegap)
    P_plasma = pfa.P(L)*1000 # in mPa

    return P_drude, P_plasma


if __name__ == "__main__":
    from sys import argv

    R = 151.3e-6
    T = 295

    omegap = 9/hbar_eV    # plasma frequency 9eV
    gamma  = 0.03/hbar_eV # dissipation 0.03
    filename = "../../materials/GoldEpsIm.dat"

    start = float(argv[1])
    stop  = float(argv[2])
    N     = int(argv[3])

    print("# L (nm), P (mPa, Drude), P (mPa, Plasma)")
    for Li in np.linspace(start, stop, N):
        pfa = PFA(R,T,get_epsilonm1(filename, omegap, gamma))
        P_drude = pfa.P(Li*1e-9)*1000 # in mPa

        pfa = PFA(R,T,get_epsilonm1(filename, omegap, gamma), omegap=omegap)
        P_plasma = pfa.P(Li*1e-9)*1000 # in mPa

        print("%.15g, %.12g, %.12g" % (Li,P_drude, P_plasma))
