import numpy as np
from math import exp, fsum, pi, log1p, isinf
from scipy import interpolate
from scipy.integrate import quad
from spence import Li2

# constants
c       = 299792458       # speed of light [m/s]
kB      = 1.38064852e-23  # Boltzmann constant [m² kg / (s² K)]
hbar    = 1.0545718e-34   # hbar [J s]
hbar_eV = 6.582119514e-16 # hbar [eV s / rad]
inf     = float("inf")    # infinity
hbarc   = hbar*c          # hbar


def __theta_plasma(a):
    """ a = c/(omegap*L)"""
    # coefficients λ_ij, Table II, [Teo, PRD 88, 045019 (2013)]
    λ = {
        (0,0): -20/pi**2+1/3,
        (1,0): +56/(3*pi**2)-32/45,
        (0,1): +56/(3*pi**2)-14/45,
        (2,0): -398/21/pi**2+401/315,
        (1,1): -796/21/pi**2+454/315,
        (0,2): -398/21/pi**2+113/315,
        (3,0): +410/21/pi**2-37/18+286/6615*pi**2,
        (2,1): +410/7/pi**2-26/7,
        (1,2): +410/7/pi**2-16/7,
        (0,3): +410/21/pi**2-79/126+pi**2/6615,
        (4,0): -69824/3465/pi**2+35141/10395-28022/99225*pi**2,
        (3,1): -279296/3465/pi**2+84176/10395-2774/14175*pi**2,
        (2,2): -139648/1155/pi**2+742/99+32/11025*pi**2,
        (1,3): -279296/3465/pi**2+43856/10395-46558/1091475*pi**2,
        (0,4): -69824/3465/pi**2+14981/10395-11962/1091475*pi**2,
        (5,0): +26732/1287/pi**2-150368/27027+4937399/5675670*pi**2-1142/63063*pi**4,
        (4,1): +133660/1287/pi**2-35026/2079+773884/567567*pi**2,
        (3,2): +267320/1287/pi**2-548024/27027+26212/51597*pi**2,
        (2,3): +267320/1287/pi**2-415724/27027+16826/81081*pi**2,
        (1,4): 133660/1287/pi**2-256888/27027+19984/81081*pi**2,
        (0,5): 26732/1287/pi**2-84218/27027+3329/62370*pi**2+8059/2522520*pi**4,
        (5,0): 26732/1287/pi**2-150368/27027+4937399/5675670*pi**2-1142/63063*pi**4
    }

    if a < 5:
        raise ValueError("a must be a>=5")

    terms = []
    maximum = 5
    for i in range(maximum+1):
        for j in range(maximum+1):
            terms.append(λ.get((i,j),0)*a**(i+j))

    s = fsum(terms)

    R = 100e-6 # exact value is irrelevant
    L = 200e-9 # exact value is irrelevant
    E_PR = -pi**3*hbarc*R/(720*L**2)

    omegap = c/(a*L)
    epsm1 = lambda xi: (omegap/xi)**2
    pfa = PFA(R, 0, epsm1)
    E = pfa.E(L)
    theta0 = E/E_PR

    theta1 = fsum(terms)/theta0
    return theta0,theta1


def theta_plasma(a):
    """Compute correction to PFA for plasma model at T=0
    E ≈ -π³ħcR/(720L²)*(θ_1 + θ_2·x), where x=L/R.
    Note that θ_1, θ_2 depend only on a.

    Parameters: a=c/(ω_P*L) (float, list, tuple or numpy array)
    Returns: θ_1, θ_2

    This function only works if a ≲ 10. For a < 5 an exception will be raised.

    Reference: Teo, PRD 88, 045019 (2013)
    """
    if type(a) == float or type(a) == int:
        return __theta_plasma(a)

    theta0 = np.zeros(len(a))
    theta1 = np.zeros(len(a))

    for j,a_j in enumerate(a):
        theta0[j], theta1[j] = __theta_plasma(a_j)

    return theta0,theta1


def psum(f, L, T, epsrel=1e-9):
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
        I,err = quad(integrand, 0, inf, epsabs=0, epsrel=epsrel)
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


def epsilonm1_from_file(filename):
    """Get dielectric function eps(xi)-1 from file"""
    with open(filename, "r") as f:
        content = f.read()

    def extract(label):
        string = "# %s" % label
        start = content.find(string)
        stop = start+content[start:].find("\n")
        line = content[start:stop]
        label,value = line.split("=")
        return float(value[:-2]) # strip eV

    omegap_low  = extract("omegap_low")  # in eV
    gamma_low   = extract("gamma_low")   # in eV
    omegap_high = extract("omegap_high") # in eV
    gamma_high  = extract("gamma_high")  # in eV

    data = np.loadtxt(filename)

    data_xi  = data[:,0]
    data_eps = data[:,1]
    xi_min, xi_max = data_xi[0], data_xi[-1]

    f = interpolate.interp1d(data_xi, data_eps, kind="linear")

    def epsilonm1(xi):
        if xi_min <= xi <= xi_max: # interpolation
            return f(xi)-1
        else: # extrapolation
            xi_ = xi/hbar_eV # in rad/s

            if xi < xi_min:
                return omegap_low**2/(xi_*(xi_+gamma_low))
            else:
                return omegap_high**2/(xi_*(xi_+gamma_high))

    return epsilonm1


class PFA:
    def __init__(self, R, T, epsm1, model="Drude", omegap=9):
        """Initialize PFA object
        R:      radius of sphere in m
        T:      temperature in Kelvin
        epsm1:  dielectric function minus 1, epsm1(xi) = eps-1
        model:  either Drude, plasma or PR (perfect reflectors)
        omegap: plasma frequency for the plasma model in eV. The value is used
                for TE Fresnel coefficient at xi=0.
        """
        self.R = R
        self.T = T
        self.epsm1 = epsm1
        self.omegap = omegap
        self.model = model


    def rp(self, xi, kappa):
        """Fresnel coefficients
        Compute Fresnel coefficients r_TE and r_TM.
        
        If epsm1(xi) is inf, r_TE=-1 and r_TM=1 will be returned for xi >= 0.
        Otherwise, if xi=0, r_TE=0 and r_TM=1 will be returned (Drude).

        Returns r_TE,r_TM
        """
        if xi == 0:
            model = self.model.lower()
            if model == "pr":
                return -1,1
            elif model == "plasma":
                omegap = self.omegap/hbar_eV
                beta = np.sqrt(1 + (omegap/(c*kappa))**2)
                rTE = (1-beta)/(1+beta)
                return rTE,1
            else: # drude
                return 0,1


        epsm1 = self.epsm1(xi)

        beta = np.sqrt(1+ (xi/(c*kappa))**2*epsm1 )

        rTE = (1-beta)/(1+beta)
        rTM = (epsm1+1-beta)/(epsm1+1+beta)

        return rTE, rTM


    def dF_xi(self, L, xi, epsabs=0, epsrel=1e-12):
        def f(t):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            rTE2, rTM2 = rTE**2, rTM**2
            exp_f = exp(-2*t)
            return t**2*exp_f*( rTE2/(1-rTE2*exp_f) + rTM2/(1-rTM2*exp_f) )

        I,err = quad(f, xi*L/c, inf, epsabs=epsabs, epsrel=epsrel)
        return I


    def dF(self, L):
        """ dF = dF/dL = - d²F/dL²"""
        return 2*self.R/L**3*psum(self.dF_xi, L, self.T)


    def F_xi(self, L, xi, epsabs=0, epsrel=1e-12):
        def f(t):
            """Integrand \sum_p t² rp² e^(-2t)/(1-rp² e^(-2t))"""
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            exp_f = exp(-2*t)
            return t*(log1p(-rTE**2*exp_f)+log1p(-rTM**2*exp_f))

        I,err = quad(f, xi*L/c, inf, epsabs=epsabs, epsrel=epsrel)
        return I


    def F(self, L):
        """Casimir force, F = -dE/dL"""
        return self.R/L**2*psum(self.F_xi, L, self.T)


    def E_xi(self, L, xi, epsabs=0, epsrel=1e-12):
        def f(t):
            kappa = t/L
            rTE, rTM = self.rp(xi,kappa)
            exp_f = exp(-2*t)
            return Li2(rTE**2*exp_f)+Li2(rTM**2*exp_f)

        I,err = quad(f, xi*L/c, inf, epsabs=epsabs, epsrel=epsrel)
        return I

    def E(self,L):
        """Casimir free energy, E"""
        return -0.5*self.R/L*psum(self.E_xi, L, self.T)


class PP:
    def __init__(self, T, epsm1, model="Drude", omegap=9):
        """Initialize PP object
        T:     temperature in Kelvin
        epsm1: dielectric function minus 1, epsm1(xi) = eps-1
        """
        self.T = T
        self.pfa = PFA(1, T, epsm1, model=model, omegap=omegap)

    def F(self, L):
        """Return F/A in N/m²"""
        return -psum(self.pfa.dF_xi, L, self.T)/(pi*L**3)

    def E(self, L):
        """Return E/A in J/m²"""
        return psum(self.pfa.F_xi, L, self.T)/(2*pi*L**2)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compute PFA in the plane-sphere geometry")
    parser.add_argument("--omegap", type=float, default=inf, help="plasma frequency in eV")
    parser.add_argument("--gamma", type=float, default=0, help="relaxation frequency in eV")
    parser.add_argument("-L", type=float, action="append", required=True, help="separation between sphere and plate in m")
    parser.add_argument("-R", type=float, required=True, help="radius of sphere in m")
    parser.add_argument("-T", type=float, default=0, help="temperature in K")
    args = parser.parse_args()

    R = args.R
    T = args.T
    omegap = args.omegap/hbar_eV # in rad/s
    gamma  = args.gamma/hbar_eV  # in rad/s

    if isinf(omegap):
        model = "PR"
        epsm1 = lambda xi: inf
    elif gamma > 0:
        model = "Drude"
        epsm1 = lambda xi: omegap**2/(xi*(xi+gamma))
    else:
        model = "Plasma"
        epsm1 = lambda xi: omegap**2/xi**2

    pfa = PFA(R, T, epsm1, model=model)

    print("# L/R, L [m], R [m], T [K], omegap [eV], gamma [eV], E_PFA*(L+R)/(hbar*c)")
    for L in args.L:
        E = pfa.E(L)
        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g" % (L/R, L, R, T, args.omegap, args.gamma, E*(L+R)/(hbar*c)))
