#!/usr/bin/python3

from math import log, pi
from glob import glob

from pyx import color, graph, text
import numpy as np
from scipy.optimize import curve_fit


class FitL:
    '''fit to c0+c1*log(x)+c2*log(x)**2+...

    deg is the highest power in log(x), so that the number of
    coefficients is given by deg+1

    '''
    def __init__(self, xvals, yvals, deg):
        self.xvals = xvals
        self.yvals = yvals
        self.deg = deg
        self.title = r'$\beta_l(%i)$' % deg

    def func(self, x, *coeffs):
        f = coeffs[0]
        for powm1, c in enumerate(coeffs[1:]):
            f = f+c*np.log(x)**(powm1+1)
        return f

    def fit(self):
        p0 = np.ones(self.deg+1)
        self.popt, pcov = curve_fit(self.func, self.xvals, self.yvals, p0)

    def expr(self):
        self.fit()
        e = 'y(x) = %.15g' % self.popt[0]
        for powm1, c in enumerate(self.popt[1:]):
            if c != 0:
                e = e + '%+.15g*x**%i' % (c, powm1+1)
        return e


class FitX:
    '''fit to c0+c1*x+c2*x**2+...

    deg is the highest power in x, so that the number of
    coefficients is given by deg+1

    '''
    def __init__(self, xvals, yvals, deg):
        self.xvals = xvals
        self.yvals = yvals
        self.deg = deg
        self.title = r'$\beta_x(%i)$' % deg

    def func(self, x, *coeffs):
        f = coeffs[0]
        for powm1, c in enumerate(coeffs[1:]):
            f = f+c*x**(powm1+1)
        return f

    def fit(self):
        p0 = np.ones(self.deg+1)
        self.popt, pcov = curve_fit(self.func, self.xvals, self.yvals, p0)

    def expr(self):
        self.fit()
        e = 'y(x) = %.15g' % self.popt[0]
        for powm1, c in enumerate(self.popt[1:]):
            if c != 0:
                e = e + '%+.15g*exp(%i*x)' % (c, powm1+1)
        return e


class FitM:
    '''fit to c0+c1*log(x)+c2*x+c3*x**2+...

    deg is the highest power, so that the number of
    coefficients is given by deg+2

    '''
    def __init__(self, xvals, yvals, deg):
        self.xvals = xvals
        self.yvals = yvals
        self.deg = deg
        self.title = r'$\beta_m(%i)$' % deg

    def func(self, x, *coeffs):
        f = coeffs[0]+coeffs[1]*np.log(x)
        for powm1, c in enumerate(coeffs[2:]):
            f = f+c*x**(powm1+1)
        return f

    def fit(self):
        p0 = np.ones(self.deg+2)
        self.popt, pcov = curve_fit(self.func, self.xvals, self.yvals, p0)

    def expr(self):
        self.fit()
        e = 'y(x) = %.15g%+.15g*x' % tuple(self.popt[0:2])
        for powm1, c in enumerate(self.popt[2:]):
            if c != 0:
                e = e + '%+.15g*exp(%i*x)' % (c, powm1+1)
        return e


class FitB:
    '''fit to c0+c1*x*log(x)+c2*x**2+c3*x**3+...

    deg is the highest power, so that the number of
    coefficients is given by deg+1

    '''
    def __init__(self, xvals, yvals, deg):
        self.xvals = xvals
        self.yvals = yvals
        self.deg = deg
        self.title = r'$\beta_b(%i)$' % deg

    def func(self, x, *coeffs):
        f = coeffs[0]+coeffs[1]*x*np.log(x)
        for powm1, c in enumerate(coeffs[2:]):
            f = f+c*x**(powm1+2)
        return f

    def fit(self):
        p0 = np.ones(self.deg+1)
        self.popt, pcov = curve_fit(self.func, self.xvals, self.yvals, p0)

    def expr(self):
        self.fit()
        e = 'y(x) = %.15g%+.15g*x*exp(x)' % tuple(self.popt[0:2])
        for powm1, c in enumerate(self.popt[2:]):
            if c != 0:
                e = e + '%+.15g*exp(%i*x)' % (c, powm1+2)
        return e


def read_data(globpattern):
    xvals = []
    yvals = []
    for filename in glob(globpattern):
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                empty = line == ""
                comment = line.startswith("#")
                if not(empty or comment):
                    LbyR, _, _, _, F = map(float, line.split(","))
                    ratio = F/pfa(LbyR)
                    beta = (ratio-1)/LbyR
                    xvals.append(LbyR)
                    yvals.append(beta)
    return np.array(xvals), np.array(yvals)


def pfa(x):
    return -pi**3/720.*(x+1)/x**2


allxvals, allyvals = read_data('eta8/*.out')
logallxvals = np.log(allxvals)
fitmin = -5
fitmax = -4.5
tobefitted = np.logical_and(fitmin <= logallxvals, logallxvals <= fitmax)
xvals = allxvals[tobefitted]
yvals = allyvals[tobefitted]
data_used = [(x, y) for x, y in zip(np.log(xvals), yvals)]
data_unused = [(x, y) for x, y in
               zip(np.log(allxvals[np.logical_not(tobefitted)]),
                   allyvals[np.logical_not(tobefitted)])]

funcdata = []
for fitclass, deg in ((FitL, 3), (FitX, 3), (FitM, 2), (FitB, 3)):
    f = fitclass(xvals, yvals, deg)
    funcdata.append(graph.data.function(f.expr(), title=f.title))
    print(fitclass.__name__)
    print(f.popt)
    print()

# linke und rechte grenze des plots
lnx_min, lnx_max = -6.2, -4
text.set(text.LatexRunner)

# do plot
ytitle = r'$\beta(x) = \left(\mathcal{F}/\mathcal{F}^\mathrm{PFA}-1\right)/x$'
g = graph.graphxy(width=12,
                  key=graph.key.key(pos="tl"),
                  x=graph.axis.lin(min=lnx_min, max=lnx_max, title=r"$x=L/R$"),
                  y=graph.axis.lin(max=-1.4, title=ytitle)
                  )
attrs = [color.gradient.RedBlue]
g.plot(funcdata, [graph.style.line()])

g.plot([
    graph.data.points(data_used,   x=1, y=2, title="used for fit"),
    graph.data.points(data_unused, x=1, y=2, title="not used for fit")],
    [graph.style.symbol(graph.style.symbol.circle, size=0.06,
                        symbolattrs=attrs)])

g.writePDFfile()
