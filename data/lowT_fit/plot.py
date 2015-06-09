#!/usr/bin/python

from __future__ import division

from math import *
from pyx import *
from numpy import polyfit
from glob import glob
from sys import exit

fit_plots = False

text.set(mode="latex")

data = []
print "# Q, L/R, F(T=0), S/(4T^3)"
for filename in glob("raw/slurm*.out"):
    fh = open(filename, "r")
    fit_x = []
    fit_x2 = []
    fit_y = []
    for line in fh:
        line = line.strip()
        if line == "" or line[0] == "#":
            continue
        
        # R/(L+R), T, F, lmax, nmax, time
        Q,T,F,lmax,nmax,time = map(float, line.split(","))
        LbyR = (1-Q)/Q

        fit_x.append(T**4)
        fit_x2.append(T)
        fit_y.append(F)
    fh.close()

    b,a = polyfit(fit_x, fit_y, 1)
    print "%.12g, %.12g, %.12g, %.12g" % (Q, LbyR, a, -b)

    if fit_plots:
        g = graph.graphxy(
            width = 8,
            x     = graph.axis.lin(title=r"$T$"),
            y     = graph.axis.lin(title=r"$\mathcal{F}$")
        )
        f = "y(x) = %.10g%+.10g*x**4" % (a,b)
        g.plot(graph.data.function(f, points=100))
        g.plot(graph.data.values(x=fit_x2, y=fit_y),
            [graph.style.symbol(graph.style.symbol.circle, size=0.1)]
        )
        s = "fits/fit_%.6f" % LbyR
        s = s.replace(".", "_")
        g.writePDFfile(s)

    data.append((LbyR,-a,-b))

data.sort()

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="tr"),
    x     = graph.axis.lin(title=r"$L/R$", max=10),
    y     = graph.axis.log(title=r"$-\mathcal{F}(T=0)$", min=1e-4)
)


g.plot(
    graph.data.points(data, x=1, y=2, title=r"$\mathcal{F}(T=0)$"),
    [graph.style.line([color.gradient.RedBlue])]
)

f_PFA = "y(x)=%.10g*(1+2*x)/(x**2+x**3)" % (pi**3/720)
g.plot(graph.data.function(f_PFA, points=1000, title=r"$\mathcal{F}_\mathrm{PFA}(T=0)$"), [graph.style.line([style.linestyle.dashed])])

f_LD = "y(x) = %.10g/(1+x)**3" % (9/(16*pi))
g.plot(graph.data.function(f_LD, points=1000, title=r"$\mathcal{F}_\mathrm{LD}(T=0)$"), [graph.style.line([style.linestyle.dashdotted])])

g.writePDFfile("F0")

g = graph.graphxy(
    width = 8,
    key = graph.key.key(pos="tl"),
    x   = graph.axis.log(title=r"$L/R$", max=10, min=0.03),
    y2  = graph.axis.lin(title=r"$-b = S/(4T^3)$", min=-7e-4, max=2e-5)
)
g.plot(
    graph.data.points(data, x=1, y=3, title=r"$S/(4T^3)$"),
    [graph.style.line([color.gradient.RedBlue])]
)

f_LD = "y(x) = %.10g/(1+x)**3" % (-9/(16*pi)/135)
g.plot(graph.data.function(f_LD, points=1000, title=r"$S_\mathrm{LD}/(4T^3)$"), [graph.style.line([style.linestyle.dashed])])

g.writePDFfile("S0")
