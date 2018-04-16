import teo
from sys import path,argv
path.append("../../../../src/python/")
from PFA import PFA, hbar, c, hbar_eV

def usage(self="data.py", exitcode=None):
    print("Usage: %s [FILES]" % self)
    print("Example: python %s data_a10/slurm-*.out" % self)
    if exitcode != None:
        exit(exitcode)


if __name__ == "__main__":
    if len(argv) <= 1:
        usage(argv[0], exitcode=1)

    files = argv[1:]

    # read data
    data = []
    for fname in files:
        with open(fname) as f:
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
    omegap = omegap_eV/hbar_eV
    a = c/(omegap*L)

    # compute prefactor of linear correction, see Teo, PRD 88, 045019 (2013)
    theta1 = teo.theta(a)

    print("# L = %.14g [m]" % L)
    print("# T = %.14g [K]" % T)
    print("# omega_p = %.14g [eV]" % omegap_eV)
    print("# omega_p*L/c = %g" % (1/a))
    print("# E = E_PFA (1+theta_1*L/R)")
    print("# theta1 = %.14g" % theta1)
    print("#")

    print("# L/R, R, ldim, E*(L+R)/(hbar*c), E_PFA*(L+R)/(hbar*c), E/E_PFA, 1-E/E_PFA, 1-E/E_PFA+theta1*L/R")

    for i in range(len(data)):
        LbyR = data[i][0]
        R = data[i][2]
        E = data[i][6]
        ldim = data[i][5]

        pfa = PFA(R, T, lambda xi: omegap**2/xi**2)
        E_PFA = pfa.E(L)*(L+R)/(hbar*c)

        print("%.12g, %.12g, %d, %.12g, %.12g, %.12g, %.12g, %.12g" % (LbyR, R, ldim, E, E_PFA, E/E_PFA, 1-E/E_PFA, E/E_PFA-1-theta1*LbyR))
