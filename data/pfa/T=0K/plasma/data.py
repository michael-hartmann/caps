from math import pi
import argparse
from sys import path
path.append("../../../../src/python/")
import PFA

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


def teo_theta1(a):
    """Compute linear correction to PFA in the plasma model
    E ≈ -π³ħcR/(720L²)*(1 + θ_1·x), x = L/R
    Parameter: a=c/(ω_P*L)

    References: Teo, PRD 88, 045019 (2013)

    Returns: θ_1
    """
    # exact values irrelevant
    R = 100e-6
    L = 250e-9

    # plasma frequency in rad/s
    omegap = PFA.c/(a*L)

    # PFA value for perfect reflectors
    E_PR = pi**3*PFA.hbarc*R/(720*L**2)

    epsm1 = lambda xi: (omegap/xi)**2
    pfa = PFA.PFA(R, 0, epsm1)
    E = pfa.E(L)
    s1 = -E/E_PR

    s2 = 0
    maximum = 5
    for i in range(maximum+1):
        for j in range(maximum+1):
            s2 += λ.get((i,j),0)*a**(i+j)

    return s2/s1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process raw data for plasma model at T=0")
    parser.add_argument("filenames", nargs="+")
    args = parser.parse_args()

    # read data
    data = []
    for filename in args.filenames:
        with open(filename) as f:
            omegap = float("nan")
            for line in f:
                line = line.strip()
                if "# omegap" in line:
                    omegap = float(line.split("=")[1])
                if line and line[0] != "#":
                    LbyR, L, R, T, ldim, E = map(float, line.split(","))
                    data.append((LbyR, L, R, T, omegap, ldim, E))

    # sort data descending
    data = list(sorted(data, reverse=True))

    # the data has constant L,T,omegap
    L = data[0][1]
    T = data[0][3]
    omegap_eV = data[0][4]

    # plasma frequency in rad/s
    omegap = omegap_eV/PFA.hbar_eV
    a = PFA.c/(omegap*L)

    # compute prefactor of linear correction, see Teo, PRD 88, 045019 (2013)
    theta1 = teo_theta1(a)

    print("# L = %.14g [m]" % L)
    print("# T = %.14g [K]" % T)
    print("# omega_p = %.14g [eV]" % omegap_eV)
    print("# omega_p*L/c = %g" % (1/a))
    print("# E = E_PFA (1+theta_1*L/R)")
    print("# theta1 = %.14g" % theta1)
    print("#")

    print("# L/R, R, ldim, E*(L+R)/(hbar*c), E_PFA*(L+R)/(hbar*c), E/E_PFA, 1-E/E_PFA, E/E_PFA-1-theta1*L/R")

    for i in range(len(data)):
        LbyR = data[i][0]
        R = data[i][2]
        E = data[i][6]
        ldim = data[i][5]

        pfa = PFA.PFA(R, T, lambda xi: omegap**2/xi**2)
        E_PFA = pfa.E(L)*(L+R)/PFA.hbarc

        print("%.12g, %.12g, %d, %.12g, %.12g, %.12g, %.12g, %.12g" % (LbyR, R, ldim, E, E_PFA, E/E_PFA, 1-E/E_PFA, E/E_PFA-1-theta1*LbyR))
