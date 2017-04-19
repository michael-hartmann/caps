import numpy as np
from math import fsum

def get_spacing(x, epsrel=1e-12):
    npts = len(x)

    # spacing of grid
    npts = len(x)
    h = (x[-1]-x[0])/(npts-1)

    # check if grid is aequidistand and ascending
    for i in range(npts-1):
        dx = x[i+1]-x[i]
        if dx < 0 or abs(1-dx/h) > epsrel:
            raise BaseException("grid must be aequidistant and ascending")

    return h,npts


def deriv2_central(y,h,i,step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    6 with stepsize step.
    """

    npts = len(y)

    if i-3*step < 0 or i+3*step >= npts:
        raise IndexError("Out of bounds")

    c = (1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90)

    terms = (
        c[0]*y[i-3*step],
        c[1]*y[i-2*step],
        c[2]*y[i-step],
        c[3]*y[i],
        c[4]*y[i+step],
        c[5]*y[i+2*step],
        c[6]*y[i+3*step]
    )

    return fsum(terms)/(step*h)**2


def deriv2_forward(y,h,i,step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    6 with stepsize step.
    """

    npts = len(y)

    if i+6*step >= npts:
        raise IndexError("Out of bounds")

    c = (203/45, -87/5, 117/4, -254/9, 33/2, -27/5, 137/180)

    terms = (
        c[0]*y[i],
        c[1]*y[i+step],
        c[2]*y[i+2*step],
        c[3]*y[i+3*step],
        c[4]*y[i+4*step],
        c[5]*y[i+5*step],
        c[6]*y[i+6*step]
    )

    return fsum(terms)/(step*h)**2


def deriv2_backward(y,h,i,step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    4 with stepsize step.
    """

    npts = len(y)

    if i-5*step < 0:
        raise IndexError("Out of bounds")

    c = (15/4, -77/6, +107/6, -13, 61/12, -5/6)

    terms = (
        c[5]*y[i-5*step],
        c[4]*y[i-4*step],
        c[3]*y[i-3*step],
        c[2]*y[i-2*step],
        c[1]*y[i-step],
        c[0]*y[i],
    )

    return fsum(terms)/(step*h)**2


def deriv2(x,y, maxstep=10, epsrel=1e-12):
    """order 4"""

    assert len(x) == len(y)

    # spacing of grid
    h = get_spacing(x, epsrel)

    # check if grid is aequidistand and ascending
    for i in range(npts-1):
        dx = x[i+1]-x[i]
        if dx < 0 or abs(1-dx/h) > epsrel:
            raise BaseException("grid must be aequidistant and ascending")


    def _d2f(f,i):
        Dk = []
        for j in range(maxstep):
            s = maxstep-j
            d = f(y,h,i,step=s)
            Dk.append( d )

        print(x[i], Dk)
        for j in range(2,maxstep):
                if abs(Dk[-3]-Dk[-2]) > abs(Dk[-2]-Dk[-1]):
                    print(s)
                    return d
        return None

    dx  = []
    dy = []

#   for i in range(npts-2*maxstep,npts):
#       d2f = _d2f(deriv2_backward,i)
#       if d2f != None:
#           dx.append(x[i])
#           dy.append(d2f)

#   for i in range(2*maxstep):
#       d2f = _d2f(deriv2_forward,i)
#       if d2f != None:
#           dx.append(x[i])
#           dy.append(d2f)

    for i in range(2*maxstep,npts-2*maxstep):
        d2f = deriv2_central(y,h,i,maxstep)
        dx.append(x[i])
        dy.append(d2f)
        
 
    return np.array(dx), np.array(dy)
