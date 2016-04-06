#!/usr/bin/python

from __future__ import division
from pyx import *
from sys import exit

N=1000
omegap = 1.32e16
gamma  = 6.9e13

drude      = "y(x) = 1 + %f/(x*(%f+x))" % (omegap**2, gamma)
drude_div  = "y(x) = %f/(x*%f)" % (omegap**2,gamma)
plasma = "y(x) = 1 + %f/(x**2)"     % (omegap**2)
plasma_div  = "y(x) = %f/(x**2)" % (omegap**2)

text.set(mode="latex")

attrs = [ graph.style.line([ style.linewidth.thick, color.gradient.RedBlue, attr.changelist([style.linestyle.solid, style.linestyle.dashed])  ]) ]

#regular = graph.axis.painter.regular(labelattrs=None)
#g = graph.graphxy(width=8,
#    x=graph.axis.lin(title=r"$\xi$ in s$^{-1}$", min=5e13, max=1.1e15),
#    y=graph.axis.lin(title=r"$\epsilon(i\xi)$",   min=0, max=1e4),
#    key=graph.key.key(pos="tr", dist=0.1)
#)

#g.plot([graph.data.function(drude,  title=r"Drude", points=N),
#        graph.data.function(plasma, title=r"Plasma", points=N)],attrs)
#xpos,  ypos = g.pos(7e13,540)
#g.text(xpos,ypos,"a)")

#g.writePDFfile("drude_lin.pdf")

px = graph.axis.painter.regular(titlepos=0.98, titledirection=None)
py = graph.axis.painter.regular(titledist=-0.3, titlepos=0.5, titledirection=None)

g = graph.graphxy(width=8,
    x=graph.axis.log(title=r"$\xi$ in s$^{-1}$", min=1e11, max=1e17),
    y2=graph.axis.log(title=r"$\epsilon(i\xi)$", min=1, max=5e7),
    key=graph.key.key(pos="tr", dist=0.1)
)


g.plot([graph.data.function(drude,       title=r"Drude", points=N),
        graph.data.function(drude_div,   title=r"$\omega_P^2/\gamma\xi$", points=N),
        graph.data.function(plasma,      title=r"Plasma", points=N),
        graph.data.function(plasma_div,  title=r"$\omega_P^2/\xi^2$", points=N)],attrs)


xpos,  ypos = g.pos(1.5e11,3)
#g.text(xpos,ypos,"b)")

g.writePDFfile("drude_log.pdf")
