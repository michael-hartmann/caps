#!/usr/bin/python3

from math import *
from pyx import *
from glob import glob
import numpy as np
from fit import *

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


polyfit = PolynomialFit(xvals, yvals)
polyfit.fitnormal(3)
f_x = polyfit.fitfuncnormal
polyfit.fitlog(3)
f_l = polyfit.fitfunclog


fit_x = xvals
fit_lnx = np.log(xvals)
fit_y = yvals

# betam
# startpunkte fuer c0,c1,c2 und c3 finden idee: wir nehmen die ersten vier
# punkte in fit_y_X und bestimmen die koefficienzen c0,c1,c2 und c3, so dass
# die fit-funktion durch diese punkte laeuft. Das benutzen wir dann als
# initialwerte fuer den fit.
M = np.zeros((4,4))
y = np.zeros(4)
for i in range(4):
    x      = fit_x[i]
    lnx    = fit_lnx[i]
    y[i]   = fit_y[i]
    M[i][0] = 1
    M[i][1] = lnx
    M[i][2] = x
    M[i][3] = x**2

c_m = np.linalg.solve(M,y)

c0 = Parameter(c_m[0])
c1 = Parameter(c_m[1])
c2 = Parameter(c_m[2])
c3 = Parameter(c_m[3])

def betam(x):
    return c0() + c1()*x + c2()*exp(x) + c3()*exp(x)**2

fit(betam, [c0, c1, c2, c3], np.array(fit_y), x=np.array(fit_lnx))
f_m = "y(x) = %.15g %+.15g*x %+.15g*exp(x) %+.15g*exp(x)**2" % (c0(), c1(), c2(), c3())


# betab
M = np.zeros((4,4))
y = np.zeros(4)
for i in range(4):
    x      = fit_x[i]
    lnx    = fit_lnx[i]
    y[i]   = fit_y[i]
    M[i][0] = 1
    M[i][1] = x*lnx
    M[i][2] = x**2
    M[i][3] = x**3

c_b = np.linalg.solve(M,y)

c0 = Parameter(c_b[0])
c1 = Parameter(c_b[1])
c2 = Parameter(c_b[2])
c3 = Parameter(c_b[3])

def betab(x):
    return c0() + c1()*exp(x)*x + c2()*exp(x)**2 + c3()*exp(x)**3

fit(betab, [c0, c1, c2, c3], np.array(fit_y), x=np.array(fit_lnx))
f_b = "y(x) = %.15g + %.15g*exp(x)*x + %.15g*exp(x)**2 + %.15g*exp(x)**3" % (c0(), c1(), c2(), c3())


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
