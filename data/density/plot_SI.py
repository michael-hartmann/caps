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

R = 150 # in 1e-6
T_min,T_max       = 0.5,2.6
LbyR_min,LbyR_max = 0.0342,1
N_LbyR,N_T        = 1000,1000

inter = Interpolate("data")
plot = []
MINIMUM = None
for i,LbyR in enumerate(linspace(LbyR_min, LbyR_max, N_LbyR)):
    for T_SI in linspace(T_min, T_max, N_T):
        LbyR, T_SI, S_SI = inter.interpolate_SI(LbyR, R*1e-6, T_SI,2)
        S_SI *= 1e26
        if S_SI > 0:
            S_SI = 0
        plot.append((LbyR, T_SI, S_SI))

        if MINIMUM == None or MINIMUM[2] > S_SI:
            MINIMUM = (LbyR, T_SI, S_SI)
    print "%.2f%%" % ((i+1)/N_LbyR*100)

print MINIMUM

z = graph.axis.linear(title=r"$S_\mathrm{SI}$ in $10^{-26}$ J/K", min=-4.5)

# R, L, T, F, lmax, nmax, time elapsed, R/L, S
g = graph.graphxy(width=10,
                  x=graph.axis.linear(title=r"$L/R$ for $R=%d\mu$m" % R),
                  y=graph.axis.linear(title=r'$T_\mathrm{SI}$ in K', min=0, max=T_max))
g.plot(graph.data.points(plot, x=1, y=2, color=3), [graph.style.density(epsilon=1e-6, gradient=gradient, coloraxis=z)])
g.dolayout()

#for i,LbyR in enumerate(linspace(LbyR_min, LbyR_max, N_LbyR)):
#    x,y = g.pos(LbyR, 539.4/(1+LbyR))
#    g.text(x,y,".")
g.writePDFfile()
