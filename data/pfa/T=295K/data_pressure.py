import numpy as np
from glob import glob
from math import pi
import casimir
from deriv import get_spacing, deriv2_central

def deriv2(x,y,step=lambda x: 1):
    """compute the second derivative. The stepsize dependent on x is given as a
    callback function"""
    d2x = []
    d2y = []

    h,npts = get_spacing(x)

    for i in range(npts):
        try:
            fpp = deriv2_central(y,h,i,step(x[i]))
            d2x.append(x[i])
            d2y.append(fpp)
        except IndexError:
            pass

    return np.array(d2x), np.array(d2y)


def stepsize(R):
    """This gives the step size dependent on separation for different radii.
    The values are hand-chosen."""

    stepsizes = {
        # radius
        151: { 300e-9: 2, 400e-9: 3, 500e-9: 4, 600e-9: 5 }
    }

    key = int(R*1e6)
    if key not in stepsizes:
        return lambda x: 1
    d = stepsizes[int(R*1e6)]

    def f(x):
        for key in reversed(sorted(d.keys())):
            if x > key:
                return d[key]
        return 1
    return f


if __name__ == "__main__":
    from sys import argv

    data = casimir.slurp(argv[1:])
    R = data[0,1]
    T = data[0,2]

    hbarc = casimir.hbar*casimir.c
    L = data[:,0]
    E_drude  = data[:,4]/(L+R)*hbarc
    E_plasma = data[:,5]/(L+R)*hbarc

    dx, dF_drude  = deriv2(L, E_drude,  step=stepsize(R))
    dx, dF_plasma = deriv2(L, E_plasma, step=stepsize(R))

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
