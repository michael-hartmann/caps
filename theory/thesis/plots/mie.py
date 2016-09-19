from __future__ import division
import numpy as np
from math import *
from scipy.special import kv,iv
from pyx import *

text.set(mode="latex")

start = 1e-4
stop = 50
N    = 500
l = 3

def a0(l):
    return pi*(-1)**l*( 2*gamma(l+1.5)-l*gamma(l+0.5) )/( l*gamma(l+0.5)**2*gamma(l+1.5) )

def b0(l):
    return pi*(-1)**(l+1)/( gamma(l+0.5)*gamma(l+1.5) )

# reflection coefficient a_l
def a(l,arg):
    return pi/2*(-1)**(l+1)* (l*iv(l+0.5,arg) - arg*iv(l-0.5,arg)) / (l*kv(l+0.5,arg) + arg*kv(l-0.5,arg))


# reflection coefficient b_l
def b(l,arg):
    return pi/2*(-1)**(l+1)* iv(l+0.5,arg)/kv(l+0.5,arg)


points = []
for x in np.logspace(log(start), log(stop), N, base=e):
    points.append((x, abs(a(l,x)), abs(b(l,x)), abs(a0(l))*(x/2)**(2*l+1), abs(b0(l))*(x/2)**(2*l+1), exp(2*x)/2, exp(2*x)/2))
    
g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="tl"),
    x     = graph.axis.log(title=r"$\chi = \frac{\xi R}{\mathrm{c}}$"),
    y     = graph.axis.log(title=r"-$a_3(i \chi)$"),
)

g.plot([
    graph.data.points(points, x=1, y=2, title=r"$a_%d(i \chi)$" % l),
    graph.data.points(points, x=1, y=4, title=r"approx. for $\chi\ll 1$"),
    graph.data.points(points, x=1, y=6, title=r"approx. for $\chi\gg 1$")],
    [graph.style.line([color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])]
)

g.writePDFfile("mie.pdf")
