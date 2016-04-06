#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from numpy import linspace

text.set(mode="latex")

def foobar(T):
    #scale = -1
    x = exp(-2*T)
    ETEE = -T**3/(4*pi)*x*(x+1)/(1-x)**3
    MTMM = -T**3/(8*pi)*x*(x+1)/(1-x)**3

    ETME = -T/(4*pi)*(T**2*x*(1+x)/(1-x)**3 + 2*T*x/(1-x)**2 + 1/(1-x) - 0.5)
    MTEM = -T/(8*pi)*(T**2*x*(1+x)/(1-x)**3 + 2*T*x/(1-x)**2 + 1/(1-x) - 0.5)

    return ETEE, MTMM, ETME, MTEM

def csch(x):
    return 1/sinh(x)

def coth(x):
    return cosh(x)/sinh(x)

def Phi(T):
    return (T*sinh(T)+cosh(T)*(T**2+sinh(T)**2))/sinh(T)**3

def dPhi(T):
    return coth(T)**2*(-(3*T**3*csch(T)**2+T))+T*(T**2+2)*csch(T)**2+coth(T)*(T**2*csch(T)**2+1)+T


f1 = []
f3 = []
for T in linspace(0.0001, 5, 1000):
    F = -3*T/(16*pi)*Phi(T)
    ETEE, MTMM, ETME, MTEM = foobar(T)
    print T,ETME
    f3.append((T,F, ETEE, MTMM, ETME, MTEM))

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="br"),
    x     = graph.axis.lin(title=r"$T$"),
    y2    = graph.axis.lin(title=r"$\mathcal{F}/\left(\frac{R}{\mathcal{L}}\right)^3$"),
)
g.plot([
    graph.data.points(f3, x=1, y=3, title=r"$E \to \mathrm{TE} \to E$"),
    graph.data.points(f3, x=1, y=4, title=r"$M \to \mathrm{TM} \to M$"),
    ],
    [graph.style.line([color.gradient.RedBlue, graph.style.line.changelinestyle])])
g.writePDFfile("tetm_change_ld")

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(pos="bl"),
    x     = graph.axis.lin(title=r"$T$"),
    y     = graph.axis.lin(title=r"$\mathcal{F}/\left(\frac{R}{\mathcal{L}}\right)^3$", max=-0.04),
)
g.plot([
    graph.data.points(f3, x=1, y=5, title=r"$E \to \mathrm{TM} \to E$"),
    graph.data.points(f3, x=1, y=6, title=r"$M \to \mathrm{TE} \to M$")
    ],
    [graph.style.line([color.gradient.RedBlue, graph.style.line.changelinestyle])])
g.writePDFfile("tetm_nochange_ld")
