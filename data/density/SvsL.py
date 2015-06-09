#!/usr/bin/python

from __future__ import division
from pyx import *
from numpy import linspace

gradient  = color.gradient.RedBlue
gradient1 = color.gradient.Jet
gradient2 = color.gradient.BlueWhite

def f(parameter):
    if parameter == 1:
        return gradient2.getcolor(parameter)
    else:
        return gradient1.getcolor(parameter)

gradient.getcolor = f

# L/R, R/(L+R), T, F, S, lmax, nmax, time
d = []
fh = open("data")
LbyR_last = None
for line in fh:
    line = line.strip()
    if line == "" or line[0] == '#':
        continue
    LbyR, RbyScriptL, T, F, S, lmax,nmax,time = map(float, line.split(","))

    if LbyR_last != LbyR:
        d.append([])

    d[-1].append((LbyR, T, S))
    LbyR_last = LbyR

def interpolate(LbyR, T):
    index_x = None
    index_y = None
    for i in range(len(d)-1):
        if d[i][0][0] <= LbyR and d[i+1][0][0] >= LbyR:
            index_x = i
            break

    for i in range(len(d[index_x])-1):
        if d[index_x][i][1] <= T and d[index_x][i+1][1] >= T:
            index_y = i
            break

    Q12 = d[index_x+0][index_y+1]
    Q22 = d[index_x+1][index_y+1]
    R2 = ( LbyR, Q12[1], Q12[2]+ (Q22[2]-Q12[2])/(Q22[0]-Q12[0])*(LbyR-Q12[0]) )

    Q11 = d[index_x+0][index_y+0]
    Q21 = d[index_x+1][index_y+0]
    R1 = ( LbyR, Q21[1], Q21[2]+ (Q11[2]-Q21[2])/(Q11[0]-Q21[0])*(LbyR-Q21[0]) )

    S = R1[2] + (R2[2]-R1[2])/(R2[1]-R1[1])*(T-R1[1])

    return [LbyR, T, S]


plot = []
for L in (0.1e-6, 0.2e-6, 0.5e-6, 1e-6, 2e-6):
for LbyR in linspace(0.045, 3.99, 1000):
    print LbyR
    #for T in linspace(0.2, 1.5, 100):
    for T_SI in linspace(60, 400, 1000):
        alpha = 0.0027438872408836173
        T = alpha*(1+LbyR)*T_SI
        p = interpolate(LbyR,T)
        p[2] *= 8.674872254535127e-23 * 1e26
        p[1] = T_SI
        if p[2] > 0:
            p[2] = 0
        plot.append(p)

z = graph.axis.linear(title=r"$S$ in $10^{-26}$ J/K", min=-4.5)

# R, L, T, F, lmax, nmax, time elapsed, R/L, S
g = graph.graphxy(width=10,
                  x=graph.axis.linear(title=r"$L/R$ for $R=1\mu$m", min=0.05, max=4),
                  y=graph.axis.linear(title=r'$T$ in K',   min=60,  max=400))
g.plot(graph.data.points(plot, x=1, y=2, color=3), [graph.style.density(epsilon=1e-6, gradient=gradient, coloraxis=z)])
g.writePDFfile()
