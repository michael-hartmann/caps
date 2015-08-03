#!/usr/bin/python

from __future__ import division

from pyx import *
from glob import glob
from math import pi,log

pfa = lambda x: -pi**3/720.*(x+1)/x**2

theta1 = 1./3 - 20/pi**2
theta2 = -4.52
f = "y(x)=1+%.15g*x+%.15g*x**2*log(x)" % (theta1,theta2)

def slurp(filename, data=[]):
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line == "" or line[0] == "#":
                continue
            # L/R, lmax, order, alpha, F(T=0)
            LbyR,lmax,order,alpha,F = map(float, line.split(","))
            data.append((LbyR, F, F/pfa(LbyR)))

    return data



data = []
for filename in glob("slurm-*.out"):
    slurp(filename, data)


text.set(mode = "latex")
g = graph.graphxy(
    width = 10,
    x = graph.axis.log(title=r"$x=L/R$", max=0.2, min=0.005),
    y = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T=0)$", min=0.8),
    key=graph.key.key(pos="tr", dist=0.1)
)

g2 = g.insert(graph.graphxy(
    width = 4,
    xpos = 1.5,
    ypos = 1,
    x = graph.axis.lin(max=0.008, min=0.005),
    y = graph.axis.lin(max=0.992, min=0.9875)
))


attrs = [color.gradient.RedBlue]
g.plot(
    graph.data.points(data, x=1, y=3, title="numerisch vs Bordag"),
    [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
)

g.plot(graph.data.function(f, title="Bimonte et al."))

g2.plot(
    graph.data.points(data, x=1, y=3, title="numerisch"),
    [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
)

g2.plot(graph.data.function(f))

g.writePDFfile("pfaT0.pdf")
