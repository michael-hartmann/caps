#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from numpy import linspace,logspace
from interpolate import Interpolate 

text.set(mode="latex")

gradient = color.gradient.WhiteBlue

def f(parameter):
    if parameter == 1:
        gradient2 = color.gradient.BlueWhite
        return gradient2.getcolor(parameter)
    else:
        gradient2 = color.gradient.ReverseJet
        return gradient2.getcolor(parameter)

gradient.getcolor = f

T_min,T_max       = 0.23,1.5
LbyR_min,LbyR_max = 0.2,5
N_LbyR,N_T        = 500,500

inter = Interpolate("data6")
plot = []
MINIMUM = None
for i,LbyR in enumerate(linspace(LbyR_min, LbyR_max, N_LbyR)):
    for T in linspace(T_min, T_max, N_T):
        LbyR, T, FbyF0 = inter.interpolate(LbyR, T, 4)
        plot.append((LbyR, T, FbyF0))

    print "%.2f%%" % ((i+1)/N_LbyR*100)


z = graph.axis.linear(title=r"$\mathcal{F}/\mathcal{F}(T=0)$")

# R, L, T, F, lmax, nmax, time elapsed, R/L, S
g = graph.graphxy(width=10,
                  x=graph.axis.linear(title=r"$L/R$"),
                  y=graph.axis.linear(title=r'$T$'))
g.plot(graph.data.points(plot, x=1, y=2, color=3), [graph.style.density(epsilon=1e-6, gradient=gradient, coloraxis=z)])
g.writePDFfile()
