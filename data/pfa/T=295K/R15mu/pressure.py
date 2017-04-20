import numpy as np
from glob import glob
from math import pi
import sys

sys.path.append("..")

import casimir
from deriv import deriv

def stepsize(R):
    if R > 550e-9:
        return 3
    if R > 450e-9:
        return 3
    if R > 300e-9:
        return 2
    return 1


if __name__ == "__main__":
    data = casimir.slurp(glob("gold_eta9.5/slurm-*.out"))
    R = data[0,1]
    T = data[0,2]

    hbarc = casimir.hbar*casimir.c
    L = data[:,0]
    E_drude  = data[:,4]/(L+R)*hbarc
    E_plasma = data[:,5]/(L+R)*hbarc

    dx, dF_drude  = deriv(L, E_drude,  deriv=2, accuracy=6, step=stepsize)
    dx, dF_plasma = deriv(L, E_plasma, deriv=2, accuracy=6, step=stepsize)

    P_drude  = 1000/(2*pi*R)*dF_drude
    P_plasma = 1000/(2*pi*R)*dF_plasma

    print("# L (m), R (m), T (K), P_Drude (mPa), P_Plasma (mPa), P_PFA_Drude (mPa), P_PFA_Plasma (mPa), "
          "P_Drude/P_PFA_Drude, P_Plasma/P_PFA_Plasma")
    for i,L in enumerate(dx):
        pfa_drude, pfa_plasma = casimir.pressure(L,R,T)
        ratio_drude  = P_drude[i]/pfa_drude
        ratio_plasma = P_plasma[i]/pfa_plasma

        print("%.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g"
              % (L, R, T, P_drude[i], P_plasma[i], pfa_drude, pfa_plasma, ratio_drude, ratio_plasma))
