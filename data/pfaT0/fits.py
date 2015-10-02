#!/usr/bin/python3

from math import *
from pyx import *
from glob import glob
import numpy as np
from fit import *

text.set(text.LatexRunner)

def pfa(x):
    return -pi**3/720.*(x+1)/x**2
    #return -pi**3/720./x**2


data_used = []
data_unused = []

# pfad zu den slurm-dateien
files = glob("eta8/*.out")

# linke und rechte grenze des plots
lnx_min, lnx_max = -6.2, -4

# bereich der punkte, die zum fitten verwendet werden
fitbereich = [-5, -4.5]
fit_x   = []
fit_lnx = []
fit_y   = []

# daten einlesen
for filename in files:
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            empty = line == ""
            comment = line.startswith("#")
            if not(empty or comment):
                # L/R, lmax, order, alpha, F(T=0)
                x,lmax,order,alpha,F = map(float, line.split(","))
                ratio = F/pfa(x)
                beta = (ratio-1)/x

                if fitbereich[0] <= log(x) <= fitbereich[1]:
                    fit_x.append(x)
                    fit_y.append(beta)
                    fit_lnx.append(log(x))
                    data_used.append((x, log(x), F, ratio, beta))
                else:
                    data_unused.append((x, log(x), F, ratio, beta))



# betal
c_l = np.polyfit(fit_lnx, fit_y, 3)
f_l = "y(x)=%.15g*x**3 + %.15g*x**2 + %.15g*x + %.15g" % (c_l[0], c_l[1], c_l[2], c_l[3])


# betax
c_x = np.polyfit(fit_x,   fit_y, 3)
f_x = "y(x)=%.15g*exp(x)**3 + %.15g*exp(x)**2 + %.15g*exp(x) + %.15g" % (c_x[0], c_x[1], c_x[2], c_x[3])


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
    graph.data.points(data_used,   x=2, y=5, title="used for fit"),
    graph.data.points(data_unused, x=2, y=5, title="not used for fit")],
    [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs)])

g.writePDFfile()

print("beta_l: ", f_l)
print("beta_x: ", f_x)
print("beta_m: ", f_m)
print("beta_b: ", f_b)
