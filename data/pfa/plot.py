#!/usr/bin/python

from __future__ import division
from pyx import *
import numpy as np
from glob import glob
from math import pi,log
import scipy.optimize as optimization

pfa = lambda x: -pi**3/720.*(x+1)/x**2

theta1 = 1./3 - 20/pi**2
theta2 = -4.52
filename = "data"
f = "y(x)=1+%.15g*x+%.15g*x**2*log(x)" % (theta1,theta2)

data = []
with open(filename, "r") as fh:
    for line in fh:
        line = line.strip()
        if line != "" and line[0] != "#":
            LbyR, T, F, error = map(float,line.split(","))
            data.append( (LbyR, T, F, F/pfa(LbyR)) )


text.set(mode = "latex")
g = graph.graphxy(
    width = 10,
    x = graph.axis.log(title=r"$x=L/R$", max=0.2, min=0.0057),
    y = graph.axis.lin(title=r"$\mathcal{F}/\mathcal{F}_\mathrm{PFA}(T=0)$", min=0.8),
    key=graph.key.key(pos="tr", dist=0.1)
)

g2 = g.insert(graph.graphxy(
    width = 4,
    xpos = 1.5,
    ypos = 1,
    x = graph.axis.lin(max=0.015, min=0.0055),
    y = graph.axis.lin(max=0.992, min=0.975)
))


attrs = [color.gradient.RedBlue]
g.plot(
    graph.data.points(data, x=1, y=4, title="numerisch vs Bordag"),
    [graph.style.symbol(graph.style.symbol.circle, size=0.04, symbolattrs=attrs)]
)

g.plot(graph.data.function(f, title="Bimonte et al."))

g2.plot(
    graph.data.points(data, x=1, y=4, title="numerisch"),
    [graph.style.symbol(graph.style.symbol.circle, size=0.06, symbolattrs=attrs)]
)

g2.plot(graph.data.function(f))

g.writePDFfile("pfa.pdf")
