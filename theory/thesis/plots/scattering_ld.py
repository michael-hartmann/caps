#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from numpy import linspace

text.set(mode="latex")

def csch(x):
    return 1/sinh(x)

def coth(x):
    return cosh(x)/sinh(x)

def Phi(T):
    return (T*sinh(T)+cosh(T)*(T**2+sinh(T)**2))/sinh(T)**3

def dPhi(T):
    return coth(T)**2*(-(3*T**3*csch(T)**2+T))+T*(T**2+2)*csch(T)**2+coth(T)*(T**2*csch(T)**2+1)+T

f1 = []
for T in linspace(0.001, 2, 500):
    F = dPhi(T)
    f1.append((T,F))

f2 = []
for T in linspace(0.001, 8, 500):
    F = dPhi(T)
    f2.append((T,F))

g1 = graph.graphxy(
    width=8,
    x  = graph.axis.lin(title=r"$T$"),
    y2 = graph.axis.lin(title=r"$S/\left(\frac{3}{16\pi} \left(\frac{R}{\mathcal{L}}\right)^3\right)$")
)
g1.plot(graph.data.points(f1, x=1, y=2), [graph.style.line([color.gradient.RedBlue])])
g1.stroke(g1.ygridpath(0), [style.linestyle.dashed])

g2 = g1.insert(graph.graphxy(
    width = 3.5,
    xpos  = 1,
    ypos  = 2.5,
    x     = graph.axis.lin(min=0, max=8),
    y2    = graph.axis.lin(min=-0.1, max=1.1))
)
g2.plot(graph.data.points(f2, x=1, y=2), [graph.style.line([color.gradient.RedBlue])])
g2.stroke(g2.ygridpath(0), [style.linestyle.dashed])

#g1.text(0.25,0.25,"b)")
g1.writePDFfile("scattering_S_ld")

f3 = []
for T in linspace(0.001, 3, 500):
    F = -T*Phi(T)
    f3.append((T,F))

g3 = graph.graphxy(
    width = 8,
    x     = graph.axis.lin(title=r"$T$"),
    y     = graph.axis.lin(title=r"$\mathcal{F}/\left(\frac{3}{16\pi} \left(\frac{R}{\mathcal{L}}\right)^3\right)$", max=-2.95),
)
g3.plot(graph.data.points(f3, x=1, y=2), [graph.style.line([color.gradient.RedBlue])])
g3.stroke(g3.ygridpath(-3), [style.linestyle.dashed])
#g3.text(0.25,0.25,"a)")
g3.writePDFfile("scattering_F_ld")
