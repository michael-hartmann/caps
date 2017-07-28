from pyx import *
import numpy as np
from sys import argv, stderr

# G. Bimonte, T. Emig, R. L. Jaffe, and M. Kardar, Casimir forces beyond the
# proximity approximation, EPL 97, 50001 (2012)
# https://dx.doi.org/10.1209/0295-5075/97/50001

theta1 = 1/3 - 20/np.pi**2

# E_PFA = -π³ħcR/(720L²)
# E_PFA*(L+R)/(ħc) = -π³/720*(1+x)/x² where x=L/R
pfa = lambda x: -np.pi**3/720.*(x+1)/x**2

filenames = argv[1:]
if len(filenames) == 0:
    print("Usage: %s FILENAMES" % argv[0], file=stderr)
    exit(1)

data = []
for filename in filenames:
    with open(filename, "r") as f:
        absrel = float("nan")

        for line in f:
            line = line.strip()
            empty = line == ""
            comment = line.startswith("#")
            if comment:
                index = line.find("absrel=")
                if index > 0:
                    line = line[index+7:]
                    absrel = float(line)
            elif not empty:
                # support old and new format
                try:
                    # L/R, lmax, order, alpha, F(T=0)
                    LbyR,lmax,order,alpha,F = map(float, line.split(","))
                except ValueError:
                    # L/R, L, R, T, ldim, F*(L+R)/(ħc)
                    LbyR,L,R,T,ldim,E = map(float, line.split(","))

                data.append((LbyR, ldim, E))

print("# LbyR, T [K], ldim, E*(L+R)/(hbar*c), E_PFA*(L+R)/(hbar*c), E/E_PFA, E/E_PFA-1-theta1*L/R")
for LbyR, ldim, E in sorted(data, key=lambda x: x[0], reverse=True):
    Epfa = pfa(LbyR)
    print("%.12g, 0, %d, %.12g, %.12g, %.12g, %.12g" % (LbyR, ldim, E, Epfa, E/Epfa, E/Epfa-1-theta1*LbyR))
