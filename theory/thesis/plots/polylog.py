#!/usr/bin/python

from __future__ import division
from pyx import *
from mpmath import polylog
from numpy import linspace

text.set(mode="latex")

start = 0
stop  = 1
N     = 200

data = []
for x in linspace(start, stop, N):
    data.append((x, polylog(2,x), polylog(3,x)))

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="tl", dist=0.1),
    x     = graph.axis.lin(title=r"$x$"),
    y     = graph.axis.lin(title=r"$f(x)$", max=1.7),
)

g.plot([
    graph.data.points(data, x=1, y=2, title=r"$f=\mathrm{Li}_2(x)$"),
    graph.data.points(data, x=1, y=3, title=r"$f=\mathrm{Li}_3(x)$"),
    graph.data.points(data, x=1, y=1, title=r"$f=x$")
    ],
    [graph.style.line([color.gradient.RedBlue])]
    )
g.writePDFfile("polylog")
