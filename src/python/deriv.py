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


def deriv_central(y, h, i, deriv=1, accuracy=6, step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    6 with stepsize step.
    """

    npts = len(y)

    coeffs = {
        (1,2): (-1/2, 0, 1/2),
        (1,4): (1/12, -2/3, 0, 2/3, -1/12),
        (1,6): (-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60),
        (1,8): (1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280),

        (2,2): (1, -2, 1),
        (2,4): (-1/12, 4/3, -5/2, 4/3, -1/12),
        (2,6): (1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90),
        (2,8): (-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560)
    }

    c = coeffs[(deriv,accuracy)]
    p = len(c)//2

    if i-p*step < 0 or i+p*step >= npts:
        raise IndexError("Out of bounds")

    terms = []
    for j in range(len(c)):
        terms.append( c[j]*y[i+(j-p)*step] )

    return fsum(terms)/(step*h)**deriv


def deriv_forward(y, h, i, accuracy=6, step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    6 with stepsize step.
    """

    raise NotImplemented("not properly implemented yet")

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


def deriv_backward(y, h, i, accuracy=6, step=1):
    """Compute the 2nd derivative of a function that is given by a vector y on
    an aequidistant grid of spacing h using central difference formula of order
    4 with stepsize step.
    """

    raise NotImplemented("not properly implemented yet")

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


def deriv(x,y, method="central", deriv=1, accuracy=6, step=lambda x: 1, factor=1):
    """compute the first or second derivative. The stepsize dependent on x is
    given as a callback function

    parameters:
    x: x values (must be equidistant)
    y: y values
    deriv: 1 for 1st derivative, 2 for 2nd derivative
    accuracy: accuracy of finite difference formula, may be 2,4,6,8
    step: callback function to determine step size
    factor: factor to multipy the derivative
    """
    d2x = []
    d2y = []

    h,npts = get_spacing(x)

    methods = { "forward": deriv_forward, "backward": deriv_backward, "central": deriv_central }
    _deriv = methods[method]

    for i in range(npts):
        try:
            fpp = _deriv(y,h,i, deriv=deriv, accuracy=accuracy, step=step(x[i]))
            d2x.append(x[i])
            d2y.append(fpp)
        except IndexError:
            pass

    return np.array(d2x), factor*np.array(d2y)
