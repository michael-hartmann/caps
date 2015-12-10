#!/usr/bin/python3

from math import log
from pyx import *
from glob import glob
import numpy as np
from scipy.optimize import curve_fit
from fit import *

class FitL:
    '''fit to c0+c1*log(x)+c2*log(x)**2+...

    deg is the highest power in log(x), so that the number of
    coefficients is given by deg+1

    '''
    def __init__(self, xvals, yvals, deg):
        self.xvals = xvals
        self.yvals = yvals
        self.deg = deg

    def func(self, x, *coeffs):
        f = coeffs[0]
        for powm1, c in enumerate(coeffs[1:]):
            f = f+c*log(x)**(powm1+1)
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

    def func(self, x, *coeffs):
        f = coeffs[0]+coeffs[1]*log(x)
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

    def func(self, x, *coeffs):
        f = coeffs[0]+coeffs[1]*x*log(x)
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

class PolynomialFit:
    def __init__(self, xvals, yvals):
        self.xvals = xvals
        self.yvals = yvals
        self.logxvals = np.log(xvals)

    def fitnormal(self, deg):
        '''fit to a polynomial in powers of x

        '''
        self.coeffsnormal = np.polyfit(self.xvals, self.yvals, deg)
        fitfunc = 'y(x) = %.15g' % self.coeffsnormal[-1]
        arg = 'exp(x)'
        for powm1, c in enumerate(self.coeffsnormal[-2::-1]):
            if c != 0:
                fitfunc = fitfunc + '%+.15g*%s**%i' % (c, arg, powm1+1)
        self.fitfuncnormal = fitfunc

    def fitlog(self, deg):
        '''fit to a polynomial in powers of log(x)

        '''
        self.coeffslog = np.polyfit(self.logxvals, self.yvals, deg)
        fitfunc = 'y(x) = %.15g' % self.coeffslog[-1]
        arg = 'x'
        for powm1, c in enumerate(self.coeffslog[-2::-1]):
            if c != 0:
                fitfunc = fitfunc + '%+.15g*%s**%i' % (c, arg, powm1+1)
        self.fitfunclog = fitfunc


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


f = FitL(xvals, yvals, 3)
f_l = f.expr()

f = FitX(xvals, yvals, 3)
f_x = f.expr()

f = FitM(xvals, yvals, 2)
f_m = f.expr()

f = FitB(xvals, yvals, 3)
f_b = f.expr()


# linke und rechte grenze des plots
lnx_min, lnx_max = -6.2, -4
text.set(text.LatexRunner)

# do plot
g = graph.graphxy(
    width = 12,
    key   = graph.key.key(pos="tl"),
    x     = graph.axis.lin(title=r"$x=L/R$", min=lnx_min, max=lnx_max),
    y     = graph.axis.lin(title=r"$\beta(x) = \left(\mathcal{F}/\mathcal{F}^\mathrm{PFA}-1\right)/x$", max=-1.4)
)
attrs = [color.gradient.RedBlue]
g.plot([
    graph.data.function(f_l, title=r"$\beta_l=a_0+a_1\log x+a_2 \log^2 x + a_3 \log^3 x$"),
    graph.data.function(f_x, title=r"$\beta_x=b_0+b_1 x+b_2 x^2 + b_3 x^3$"),
    graph.data.function(f_m, title=r"$\beta_m=c_0+c_1\log x + c_2 x + c_3 x^2$"),
    graph.data.function(f_b, title=r"$\beta_b=d_0+d_1 x \log x + d_2 x^2 + d_3 x^3$")
    ], [graph.style.line()]
)

g.plot([
    graph.data.points(data_used,   x=1, y=2, title="used for fit"),
    graph.data.points(data_unused, x=1, y=2, title="not used for fit")],
    [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs)])

g.writePDFfile()

print("beta_l: ", f_l)
print("beta_x: ", f_x)
print("beta_m: ", f_m)
print("beta_b: ", f_b)
