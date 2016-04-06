#!/usr/bin/python

from scipy.special import iv,kv
from pyx import *
import numpy as np

text.set(mode="latex")

g = graph.graphxy(width=8, key=graph.key.key(pos="tl", dist=0),
    x=graph.axis.linear(title="$x$"),
    y=graph.axis.linear(title=r"$I_\nu(x)$", max=5)
)
liste2 = []
attrs = [ graph.style.line([ style.linewidth.thick, color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed])  ]) ]
d = []

for l,t in ((0, r"$\nu=0$"), (0.5, r"$\nu=\frac{1}{2}$"), (1, r"$\nu=1$"), (1.5, r"$\nu=\frac{3}{2}$"), (2, r"$\nu=2$"), (2.5, r"$\nu=\frac{5}{2}$")):
    liste = []
    for x in np.linspace(0, 4, 100):
        liste.append((x,iv(l,x)))
    d.append(graph.data.points(liste, x=1, y=2, title=t))
    liste2.append(liste)

g.plot(d, attrs)
xpos, ypos = g.pos(3.7,0.3)
#g.text(xpos,ypos,"a)")

g.writePDFfile("iv.pdf")

# kv

g = graph.graphxy(width=8, key=graph.key.key(pos="tr", dist=0),
    x=graph.axis.linear(title="$x$"),
    y2=graph.axis.linear(title=r"$K_\nu(x)$", max=5)
)
liste2 = []
d = []

for l,t in ((0, r"$\nu=0$"), (0.5, r"$\nu=\frac{1}{2}$"), (1, r"$\nu=1$"), (1.5, r"$\nu=\frac{3}{2}$"), (2, r"$\nu=2$"), (2.5, r"$\nu=\frac{5}{2}$")):
    liste = []
    for x in np.linspace(0.001, 2.5, 100):
        liste.append((x,kv(l,x)))
    d.append(graph.data.points(liste, x=1, y=2, title=t))
    liste2.append(liste)

g.plot(d, attrs)
xpos, ypos = g.pos(0.1,0.3)
#g.text(xpos,ypos,"b)")

g.writePDFfile("kv.pdf")
