#!/usr/bin/python

from math import *
from scipy.special import iv,kv
from pyx import *
import numpy as np

text.set(mode="latex")

xmin = 0
xmax = 16

g = graph.graphxy(width=8, key=graph.key.key(pos="tl", dist=0),
    x=graph.axis.lin(title="$x$", min=0.1, max=xmax),
    y=graph.axis.log(title="$y$", min=1e-3, max=1e7)
)
attrs = [ graph.style.line([ style.linewidth.thick, color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])  ]) ]

f = []
l=2.5
for x in np.linspace(xmin, xmax, 100):
    f.append((x,iv(l,x)))

g.plot([
    graph.data.points(f, x=1, y=2, title=r"$y=I_\frac{3}{2}(x)$"),
    graph.data.function("y(x)=%f*(x/2)**%f" % (1/gamma(l+1),l), min=xmin, max=xmax, title=r"$y=\frac{1}{\Gamma(1+\frac{3}{2})} \left(\frac{x}{2}\right)^\frac{3}{2}$"),
    graph.data.function("y(x)=exp(x)/sqrt(2*pi*x)",       min=xmin, max=xmax, title=r"$y=\frac{e^x}{\sqrt{2\pi x}}$")
], attrs)
xpos, ypos = g.pos(15,5e-3)
#g.text(xpos,ypos,"a)")

g.writePDFfile("iv_approx.pdf")

xmin = 0.1
xmax = 16

g = graph.graphxy(width=8, key=graph.key.key(pos="tr", dist=0),
    x=graph.axis.lin(title="$x$", min=xmin, max=xmax),
    y2=graph.axis.log(title="$y$", min=1e-9)
)

f = []
l=2.5
for x in np.linspace(xmin, xmax, 100):
    f.append((x,kv(l,x)))

g.plot([
    graph.data.points(f, x=1, y=2, title=r"$y=K_\frac{3}{2}(x)$"),
    graph.data.function("y(x)=%f*(2/x)**%f" % (gamma(l)/2,l), min=xmin, max=xmax, title=r"$y=\frac{1}{2}\, \Gamma(\frac{3}{2}) \left(\frac{2}{x}\right)^\frac{3}{2}$"),
    graph.data.function("y(x)=sqrt(pi/(2*x))*exp(-x)",       min=xmin, max=xmax, title=r"$y=\sqrt{\frac{\pi}{2x}} e^{-x}$")
], attrs)
xpos, ypos = g.pos(0.7,7e-9)
#g.text(xpos,ypos,"b)")

g.writePDFfile("kv_approx.pdf")
