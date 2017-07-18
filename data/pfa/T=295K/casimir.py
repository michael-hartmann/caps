import numpy as np
from math import exp, fsum, pi, sqrt, log1p
from scipy.integrate import quad
from scipy import interpolate
from os.path import dirname
from deriv import deriv
from sys import stdout

# constants
c       = 299792458       # speed of light [m/s]
kB      = 1.38064852e-23  # Boltzmann constant [m² kg / (s² K)]
hbar    = 1.0545718e-34   # hbar [J s]
hbar_eV = 6.582119514e-16 # hbar [eV s]

def slurp(filenames):
    hbarc = hbar*c # hbar*c
    data = []

    for i,filename in enumerate(filenames):
        with open(filename, "r") as f:
            drude = plasma = np.nan
            for line in f:
                line = line.strip()
                if "# plasma" in line:
                    _,line = line.split("=", 1)

                    line = line.strip()
                    plasma = float(line[:line.find(" ")])
                    continue
                if "# xi=0" in line:
                    line = line[16:]
                    line = line[:line.find(",")]
                    drude = float(line)
                    continue
                if line == "" or line[0] == "#":
                    continue

                # L/R, L, R, T, ldim, E*(L+R)/(ħc)
                LbyR, L, R, T, ldim, E_scaled = map(float, line.split(","))

                # free energy (drude)
                E_drude = E_scaled*hbarc/(L+R)

                # difference free energy plasma-drude; i.e., E_plasma-E_drude
                # Please note that we don't have kb*T/2 here, because of the
                # summation over m which gives a factor of 2.
                delta_E_plasma = (plasma-drude)*kB*T

                # free energy (plasma)
                E_plasma = E_drude+delta_E_plasma

                # L in m, R in m, T in K, E_drude and E_plasma in J
                data.append((L,R,T,ldim,E_drude,E_plasma))

    return np.array(sorted(data))


def get_epsilonm1(filename):
    data = np.loadtxt(filename)

    data_xi  = data[:,0]
    data_eps = data[:,1]
    xi_min, xi_max = data_xi[0], data_xi[-1]

    f = interpolate.interp1d(data_xi, data_eps, kind="linear")

    def epsilonm1(xi,args):
        assert xi_min <= xi <= xi_max
        return f(xi)-1

    return epsilonm1


class PFA:
    def __init__(self, R, T, filename=None, args=None):
        """R in m, T in K"""
        if filename == None:
            filename = dirname(__file__) + "/../../materials/GoldDalvit.dat"
        self.R = R
        self.T = T
        self.epsilonm1 = get_epsilonm1(filename)
        self.args = args


    def rp(self, xi, kappa, omegap):
        if xi == 0:
            # plasma
            if omegap > 0:
                beta = sqrt(1+(omegap/(c*kappa))**2)
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


    def dF_xi(self, L, xi, omegap=0, epsabs=0, epsrel=1e-10):
        def f(t,xi):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa,omegap)
            rTE2, rTM2 = rTE**2, rTM**2
            exp_f = exp(-2*t)
            return t**2*exp_f*( rTE2/(1-rTE2*exp_f) + rTM2/(1-rTM2*exp_f) )

        I,err = quad(f, xi*L/c, np.inf, args=(xi,), epsabs=epsabs, epsrel=epsrel)
        return I

    def dF(self, L, omegap):
        """ dF = dF/dL = - d²F/dL²"""
        plasma_xi0 = self.dF_xi(L,0,omegap=omegap)

        terms = []
        n = 0
        while True:
            xi_n = 2*pi*n*kB*self.T/hbar
            I = self.dF_xi(L, xi_n)
            terms.append(I)

            if abs(I/terms[0]) < 1e-10:
                break

            n += 1

        terms[0] /= 2
        drude = (kB*self.R*self.T)/L**3*2*fsum(terms)

        terms[0] = plasma_xi0/2
        plasma = (kB*self.R*self.T)/L**3*2*fsum(terms)
        return drude,plasma


    def F_xi(self, L, xi, omegap=0, epsabs=0, epsrel=1e-10):
        def f(t,xi):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa,omegap)
            exp_f = exp(-2*t)
            return t*(log1p(-rTE**2*exp_f)+log1p(-rTM**2*exp_f))

        I,err = quad(f, xi*L/c, np.inf, args=(xi,), epsabs=epsabs, epsrel=epsrel)
        return I

    def F(self, L, omegap):
        """Casimir force, F = -dE/dL"""
        plasma_xi0 = self.F_xi(L,0,omegap=omegap)

        T = self.T
        n = 0
        terms = []
        while True:
            xi_n = 2*pi*n*kB*T/hbar
            I = self.F_xi(L, xi_n)
            terms.append(I)

            if abs(I/terms[0]) < 1e-10:
                break

            n += 1

        terms[0] /= 2
        drude = (kB*T*self.R/L**2)*fsum(terms)

        terms[0] = plasma_xi0/2
        plasma = (kB*T*self.R/L**2)*fsum(terms)
        return drude,plasma



def pfa_gradient(L,R,T):
    omegap = 9/hbar_eV # plasma frequency 9eV

    pfa = PFA(R,T)
    dF_drude, dF_plasma = pfa.dF(L,omegap) # in Pa

    return dF_drude, dF_plasma


def pfa_force(L,R,T):
    omegap = 9/hbar_eV # plasma frequency 9eV

    pfa = PFA(R,T)
    return pfa.F(L,omegap)


def print_force(files, stepsize, f=stdout):
    data = slurp(files)
    R = data[0,1]
    T = data[0,2]

    L = data[:,0]
    E_drude  = data[:,4]
    E_plasma = data[:,5]

    E_plasma_xi0 = E_plasma-E_drude

    dx, F_drude      = deriv(L, E_drude,      deriv=1, accuracy=6, step=stepsize)
    dx, F_plasma     = deriv(L, E_plasma,     deriv=1, accuracy=6, step=stepsize)
    dx, F_plasma_xi0 = deriv(L, E_plasma_xi0, deriv=1, accuracy=6, step=stepsize)

    F_drude *= -1
    F_plasma *= -1
    F_plasma_xi0 *= -1

    print("# L (m), R (m), T (K), F_Drude (N), F_Plasma (N), F_PFA_Drude (N), F_PFA_Plasma (N), "
          "F_Drude/F_PFA_Drude, F_Plasma/F_PFA_Plasma, F_Plasma_xi0 (N), F_Plasma_xi0/F_PFA_Plasma_xi0", file=f)
    for i,L in enumerate(dx):
        pfa_drude, pfa_plasma = pfa_force(L,R,T)
        ratio_drude  = F_drude[i]/pfa_drude
        ratio_plasma = F_plasma[i]/pfa_plasma

        ratio_plasma_xi0 = F_plasma_xi0[i]/(pfa_plasma-pfa_drude)

        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g"
               % (L, R, T, F_drude[i], F_plasma[i], pfa_drude, pfa_plasma, ratio_drude, ratio_plasma, F_plasma_xi0[i], ratio_plasma_xi0), file=f)



def print_gradient(files, stepsize, f=stdout):
    data = slurp(files)
    R = data[0,1]
    T = data[0,2]

    L = data[:,0]
    E_drude  = data[:,4]
    E_plasma = data[:,5]

    E_plasma_xi0 = E_plasma-E_drude

    dx, dF_drude      = deriv(L, E_drude,      deriv=2, accuracy=6, step=stepsize)
    dx, dF_plasma     = deriv(L, E_plasma,     deriv=2, accuracy=6, step=stepsize)
    dx, dF_plasma_xi0 = deriv(L, E_plasma_xi0, deriv=2, accuracy=6, step=stepsize)

    dF_drude *= -1
    dF_plasma *= -1
    dF_plasma_xi0 *= -1

    print("# L (m), R (m), T (K), F'_Drude (mPa), F'_Plasma (N/m), F'_PFA_Drude (N/m), F'_PFA_Plasma (N/m), "
          "F'_Drude/F'_PFA_Drude, F'_Plasma/F'_PFA_Plasma, F'_Plasma_xi0 (N/m), F'_Plasma_xi0/F'_PFA_Plasma_xi0", file=f)
    for i,L in enumerate(dx):
        pfa_drude, pfa_plasma = pfa_gradient(L,R,T)
        ratio_drude  = dF_drude[i]/pfa_drude
        ratio_plasma = dF_plasma[i]/pfa_plasma

        ratio_plasma_xi0 = dF_plasma_xi0[i]/(pfa_plasma-pfa_drude)

        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g"
              % (L, R, T, dF_drude[i], dF_plasma[i], pfa_drude, pfa_plasma, ratio_drude, ratio_plasma, dF_plasma_xi0[i], ratio_plasma_xi0), file=f)
