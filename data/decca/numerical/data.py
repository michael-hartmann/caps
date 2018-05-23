import numpy as np
from sys import path,argv

path.append("../../../src/python/")

from deriv import deriv
from PFA import *

def slurp(filename):
    with open(filename) as f:
        E_drude_xi0  = float("nan")
        E_plasma_xi0 = float("nan")
        for line in f:
            line = line.strip()
            if "# plasma" in line:
                s = line
                E_plasma_xi0 = float(s.split("=",1)[1].strip().split(" ", 1)[0])
            elif "# xi=0, logdetD=" in line:
                s = line
                E_drude_xi0 = float(s.split(",")[1].strip().split("=")[1])
            elif len(line) and line[0] != "#":
                LbyR, L, R, T, ldim, E_ = map(float, line.split(","))
                E_drude = E_*hbarc/(L+R)
                E_plasma = E_drude+(E_plasma_xi0-E_drude_xi0)*kB*T
                return R, L, E_drude, E_plasma


if __name__ == "__main__":
    T = 295 # temperature is 295K

    drude = []
    plasma = []
    for filename in argv[1:]:
        R, L, E_drude, E_plasma = slurp(filename)
        drude .append((L,E_drude))
        plasma.append((L,E_plasma))

    drude  = np.array(sorted(drude))
    plasma = np.array(sorted(plasma))

    L_drude,  F_drude  = deriv(drude[:,0],  drude[:,1],  accuracy=6, step=lambda x:1, factor=-1)
    L_plasma, F_plasma = deriv(plasma[:,0], plasma[:,1], accuracy=6, step=lambda x:1, factor=-1)

    epsm1 = epsilonm1_from_file("../../materials/GoldDalvit.dat")
    pfa_drude  = PFA(R,T,epsm1,model="Drude")
    pfa_plasma = PFA(R,T,epsm1,model="plasma")

    print("# R = %gnm" % (R*1e9))
    print("#")
    print("# L (nm), F_Drude/F_PR, F_PFA_Drude/F_PR, F_plasma/F_PR, F_PFA_plasma/F_PR")
    for i,Li in enumerate(L_drude):
        F_PR = -np.pi**3*hbarc*R/(360*Li**3)

        F_pfa_drude  = pfa_drude.F(Li)
        F_pfa_plasma = pfa_plasma.F(Li)
        print("%g, %g, %g, %g, %g" % (Li*1e9, F_drude[i]/F_PR, F_pfa_drude/F_PR, F_plasma[i]/F_PR, F_pfa_plasma/F_PR))
