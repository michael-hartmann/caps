import numpy as np
from glob import glob
from math import pi
import sys

sys.path.append("..")

import casimir
from deriv import deriv

def stepsize(R):
    if R > 300e-9:
        return 2
    return 1


if __name__ == "__main__":
    hbarc = casimir.hbar*casimir.c

    data = casimir.slurp(glob("gold_eta10/slurm-*.out"))
    R = data[0,1]
    T = data[0,2]

    L = data[:,0]
    E_drude  = data[:,4]/(L+R)*hbarc
    E_plasma = data[:,5]/(L+R)*hbarc

    dx, F_drude  = deriv(L, E_drude,  deriv=1, accuracy=6, step=stepsize)
    dx, F_plasma = deriv(L, E_plasma, deriv=1, accuracy=6, step=stepsize)

    F_drude *= -1
    F_plasma *= -1

    print("# L (m), R (m), T (K), F_Drude, F_Plasma, F_PFA_Drude, F_PFA_Plasma, "
          "F_Drude/F_PFA_Drude, F_Plasma/F_PFA_Plasma")
    for i,L in enumerate(dx):
        pfa_drude, pfa_plasma = casimir.force(L,R,T)
        ratio_drude  = F_drude[i]/pfa_drude
        ratio_plasma = F_plasma[i]/pfa_plasma

        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g"
               % (L, R, T, F_drude[i], F_plasma[i], pfa_drude, pfa_plasma, ratio_drude, ratio_plasma))
