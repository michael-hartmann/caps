#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *

text.set(mode="latex")

c = canvas.canvas()

attr = [ style.linewidth.Thin, color.rgb.black, deco.filled([color.cmyk.Periwinkle]) ]
x,y,r = 0,0,1
c.stroke(path.circle(x,y,r), attr)
c.stroke(path.rect(-4.0,-3.6,8,0.1), attr)

for phi in [-115, -135, -155]:
    phi = radians(phi)
    x0,y0 = 1.1*cos(phi), 1.1*sin(phi)
    x1,y1 = 1.6*cos(phi), 1.6*sin(phi)
    c.stroke(path.line(x0,y0, x1,y1), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(-2.7, -0.3, r"$|\ell_2,m,P_2\rangle$"))
c.insert(text.text(-1.5, -1.6, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}"))

for phi in [-115, -135, -155]:
    phi = radians(phi+90)
    x0,y0 = 1.1*cos(phi), 1.1*sin(phi)
    x1,y1 = 1.6*cos(phi), 1.6*sin(phi)
    c.stroke(path.line(x0,y0, x1,y1), [ deco.barrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(1.4, -0.3, r"$|\ell_1,m,P_1\rangle$"))
c.insert(text.text(1.2, -1.6, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {5}}}"))
c.insert(text.text(-0.2, -1.35, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {6}}}"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = -2, -2
    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(-3.3, -1.9, r"$|\mathbf{k},p,-\rangle$"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = -1.5, -3.3
    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(-2.8, -3.2, r"$|\mathbf{k},p,-\rangle$"))
c.insert(text.text(-1.5, -2.8, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = 1.5,-2
    c.stroke(path.line(x0+i*deltax,y0-deltay, x0+i*deltax+0.2,y0+deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(2.4, -1.9, r"$|\mathbf{k},p,+\rangle$"))
c.insert(text.text(1.3, -2.8, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = -1.5, -3.3
    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(-2.8, -3.2, r"$|\mathbf{k},p,-\rangle$"))
c.insert(text.text(-0.2, -3.4, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = 1, -3.3
    c.stroke(path.line(x0+i*deltax,y0-deltay, x0+i*deltax+0.2,y0+deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(1.8, -3.2, r"$|\mathbf{k},p,+\rangle$"))
#c.insert(text.text(1.8, -2.8, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}"))

#for i in range(3):
#    deltay = 0.25
#    deltax = 0.25
#    x0,y0 = -0.5-0.75/2, -2.65
#    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
#for i in range(3):
#    deltay = 0.25
#    deltax = 0.25
#    x0,y0 = 0.5-0.75/2, -2.65
#    c.stroke(path.line(x0+i*deltax+0.2,y0+deltay, x0+i*deltax+0.0,y0-deltay), [ deco.barrow([deco.stroked([deco.stroked.clear])])])
##c.insert(text.text(-0.25, -2.3, r"$(2)$"))
#
#for i in range(3):
#    deltay = 0.25
#    deltax = 0.25
#    x0,y0 = 2-0.4, -2
#    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.barrow([deco.stroked([deco.stroked.clear])])])
##c.insert(text.text(+2., -2.6, r"$(3)$"))
#
#for phi in [-135+90, -150+90, -120+90]:
#    phi = radians(phi)
#    x0,y0 = 1.1*cos(phi), 1.1*sin(phi)
#    x1,y1 = 1.6*cos(phi), 1.6*sin(phi)
#    c.stroke(path.line(x0,y0, x1,y1), [ deco.barrow([deco.stroked([deco.stroked.clear])])])

c.writePDFfile("roundtrip.pdf")
