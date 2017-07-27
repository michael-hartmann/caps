import numpy as np
from glob import glob
import sys

sys.path.append("../../../../src/python/")

from PFA import PFA, hbar, hbar_eV, c

data = []
filenames = sorted(glob("eta10/slurm-*.out"))
for filename in filenames:
    omegap,gamma = np.nan, np.nan
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if "# omegap" in line:
                omegap = float(line[11:])
            elif "# gamma" in line:
                gamma = float(line[10:])
            elif len(line) > 0 and line[0] != "#":
                LbyR, L, R, T, ldim, E = map(float, line.split(","))

                epsm1 = lambda xi: (omegap/hbar_eV)**2/(xi*(xi+gamma/hbar_eV))
                pfa = PFA(R,T,epsm1)
                Epfa = pfa.E(L) * (L+R)/(hbar*c)

                data.append((LbyR, T, omegap, gamma, ldim, E, Epfa, E/Epfa))


print("# LbyR, T [K], omegap [eV], gamma [eV], ldim, E*(L+R)/(hbar*c), E_PFA*(L+R)/(hbar*c), E/E_PFA")
for p in sorted(data, reverse=True):
    print("%.12g, %g, %g, %g, %d, %.12g, %.12g, %.12g" % p)
