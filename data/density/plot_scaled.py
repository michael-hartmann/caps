#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *
from numpy import linspace
from interpolate import Interpolate 

gradient = color.gradient.WhiteBlue

def f(parameter):
    if parameter == 1:
        gradient2 = color.gradient.BlueWhite
        return gradient2.getcolor(parameter)
    else:
        gradient2 = color.gradient.ReverseJet
        return gradient2.getcolor(parameter)

gradient.getcolor = f

T_min,T_max       = 0.207,1.5
LbyR_min,LbyR_max = 0.03415,1

N_LbyR,N_T        = 1000,1000

inter = Interpolate("data")
plot = []
MINIMUM = None
for i,LbyR in enumerate(linspace(LbyR_min, LbyR_max, N_LbyR)):
    for T in linspace(T_min, T_max, N_T):
        LbyR, T, S = inter.interpolate(LbyR, T, 2)
        if S > 0:
            S = 0
        plot.append((LbyR, T, S))

        if MINIMUM == None or MINIMUM[2] > S:
            MINIMUM = (LbyR, T, S)
    print "%.2f%%" % ((i+1)/N_LbyR*100)

print MINIMUM

z = graph.axis.linear(title=r"$S$")

# R, L, T, F, lmax, nmax, time elapsed, R/L, S
g = graph.graphxy(width=10,
                  x=graph.axis.linear(title=r"$L/R$", min=0, max=LbyR_max),
                  y=graph.axis.linear(title=r'$T$', min=0, max=T_max))
g.plot(graph.data.points(plot, x=1, y=2, color=3), [graph.style.density(epsilon=1e-6, gradient=gradient, coloraxis=z)])
g.writePDFfile()
