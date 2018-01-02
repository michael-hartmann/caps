import numpy as np
from glob import glob

import sys
sys.path.append("../../../../src/python/")
import PFA
import deriv

def slurp(filenames):
    data = []
    for filename in filenames:
        with open(filename, "r") as f:
            omegap = float("nan")
            gamma  = float("nan")
            for line in f:
                if len(line) == 0:
                    continue
                if "=" in line:
                    left,right = map(str.strip, line[1:].strip().split("=",1))
                    if left == "omegap":
                        omegap = float(right)
                    elif left == "gamma":
                        gamma = float(right)
                if line[0] != "#":
                    LbyR, L, R, T, ldim, E_ = map(float, line.split(","))
                    # L, R, T, F*(L+R)/(ħc), omegap, gamma
                    data.append((L,R,T,E_,omegap,gamma))

    return np.array(sorted(data))


if __name__ == "__main__":
    # constants
    hbarc = PFA.hbar*PFA.c
    hbar_eV = PFA.hbar_eV

    # read data
    data = slurp(glob("raw/slurm-*.out"))

    R = data[0,1]  # radius in m
    T = data[0,2]  # temperature in m
    omegap = data[0,4] # plasma frequency in eV
    gamma  = data[0,5] # relaxation frequency in eV

    L = data[:,0] # surface-to-surface separation in m
    E = data[:,3]*hbarc/(L+R) # free energy in J

    dL, dE   = deriv.deriv(L,E, deriv=1, accuracy=6, step=lambda x: 1) # dE/dL
    d2L, d2E = deriv.deriv(L,E, deriv=2, accuracy=6, step=lambda x: 1) # d²E/dL²

    header  = "# R = %g (µm)\n" % (R*1e6)
    header += "# T = %g (K)\n" % T
    header += "# model = Drude\n"
    header += "# omegap = %g (eV)\n" % omegap
    header += "# gamma = %g (eV)\n" % gamma

    # dielectric function minus 1, eps(xi)-1
    epsm1 = lambda xi: (omegap/hbar_eV)**2/(xi*(xi+gamma/hbar_eV))

    # pfa object
    pfa = PFA.PFA(R, T, epsm1)

    with open("energy.csv", "w") as f:
        print(header, file=f)
        print("# L (nm), E (J), E/Epfa", file=f)
        for i,L_ in enumerate(L):
            Epfa = pfa.E(L_)
            print("%.15g, %.15g, %.15g" % (L_*1e9, E[i], E[i]/Epfa), file=f) 

    with open("force.csv", "w") as f:
        print(header, file=f)
        print("# L (nm), F (N), F/Fpfa", file=f)
        for i,L_ in enumerate(dL):
            Fpfa = pfa.F(L_)
            print("%.15g, %.15g, %.15g" % (L_*1e9, -dE[i], -dE[i]/Fpfa), file=f) 

    with open("gradient.csv", "w") as f:
        print(header, file=f)
        print("# L (nm), F' (N/m), F'/F'pfa", file=f)
        for i,L_ in enumerate(d2L):
            dFpfa = pfa.dF(L_)
            print("%.15g, %.15g, %.15g" % (L_*1e9, -d2E[i], -d2E[i]/dFpfa), file=f) 
