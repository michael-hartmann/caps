import argparse
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
    parser = argparse.ArgumentParser(description="Process raw data")
    parser.add_argument("--order",    help="order of finite difference formula (2,4,6 or 8)", type=int, default=8)
    parser.add_argument("--stepsize", help="stepsize", type=int, default=1)
    parser.add_argument("filenames", nargs="+")
    args = parser.parse_args()
    order    = args.order
    stepsize = args.stepsize

    T = 295 # temperature is 295K

    drude = []
    plasma = []
    for filename in args.filenames:
        R, L, E_drude, E_plasma = slurp(filename)
        drude .append((L,E_drude))
        plasma.append((L,E_plasma))

    drude  = np.array(sorted(drude))
    plasma = np.array(sorted(plasma))

    L_drude,  dF_drude  = deriv(drude[:,0],  drude[:,1],  accuracy=order, deriv=2, step=lambda x: stepsize, factor=-1)
    L_plasma, dF_plasma = deriv(plasma[:,0], plasma[:,1], accuracy=order, deriv=2, step=lambda x: stepsize, factor=-1)

    epsm1 = epsilonm1_from_file("../../materials/GoldDalvit.dat")
    pfa_drude  = PFA(R,T,epsm1,model="Drude")
    pfa_plasma = PFA(R,T,epsm1,model="plasma")

    print("# optical data by Diego Dalvit, omegap = 9eV")
    print("# T = %gK" % T)
    print("# R = %gnm" % (R*1e9))
    print("#")
    print("# P_PR = pi³*hbar*c*R/(120*L^4)")
    print("#")
    print("# order of finite difference formula: %d" % order)
    print("# stepsize: %d" % stepsize)
    print("#")
    print("# L (nm), P_Drude/F_PR, P_PFA_Drude/P_PR, P_plasma/P_PR, P_PFA_plasma/P_PR")
    for i,Li in enumerate(L_drude):
        dF_PR = np.pi**3*hbarc*R/(120*Li**4)

        dF_pfa_drude  = pfa_drude.dF(Li)
        dF_pfa_plasma = pfa_plasma.dF(Li)
        print("%g, %.8g, %.8g, %.8g, %.8g" % (Li*1e9, dF_drude[i]/dF_PR, dF_pfa_drude/dF_PR, dF_plasma[i]/dF_PR, dF_pfa_plasma/dF_PR))
