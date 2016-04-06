#!/usr/bin/python

from __future__ import division
from math import *
from pyx import *

text.set(mode="latex")

c = canvas.canvas()

hoehe = 3.5
breite = 2.5*hoehe
dicke = 0.1

c.insert(text.text(0.05, hoehe+0.35, r"$z$"))
c.stroke(path.line(-0.1, 0, -0.1, hoehe+0.45), [ deco.earrow([deco.stroked([deco.stroked.clear])])])

attr = [ style.linewidth.Thin, color.rgb.black, deco.filled([color.cmyk.Periwinkle]) ]
attr2 = [
    deco.earrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
    deco.barrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
]

c.stroke(path.rect(0, hoehe, breite, 0.1), attr)
c.stroke(path.rect(0, 0, breite, 0.1), attr)

c.stroke(path.line(breite/14, dicke, breite/14, hoehe), attr2)
c.insert(text.text(breite/35, hoehe/2, r"$d$"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = 1,hoehe/2
    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(0.8, 0.5+hoehe/2, r"$|\mathbf{k},p,-\rangle$"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = breite-1,hoehe/2
    c.stroke(path.line(x0+i*deltax+0.2,y0-deltay, x0+i*deltax,y0+deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(0.8, 0.5+hoehe/2, r"$|\mathbf{k},p,-\rangle$"))
c.insert(text.text(1.8, hoehe/2, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}"))
c.insert(text.text(7.2, hoehe/2, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}"))


for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = breite/2-0.5, deltay+dicke
    c.stroke(path.line(x0+i*deltax,y0+deltay, x0+i*deltax+0.2,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(breite/2-1.9, 0.5, r"$|\mathbf{k},p,-\rangle$"))
c.insert(text.text(breite/2+1.4, 0.5, r"$|\mathbf{k},p,+\rangle$"))
c.insert(text.text(breite/2+0.15, 0.7, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}"))
c.insert(text.text(breite/2+0.15, hoehe-0.8, r"\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = breite/2+0.5, hoehe-deltay
    c.stroke(path.line(x0+i*deltax+0.2,y0-deltay, x0+i*deltax,y0+deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = breite/2-0.5, hoehe-deltay
    c.stroke(path.line(x0+i*deltax+0.2,y0+deltay, x0+i*deltax,y0-deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(breite/2-1.9, hoehe-0.4, r"$|\mathbf{k},p,-\rangle$"))
c.insert(text.text(breite/2+1.4,   hoehe-0.4, r"$|\mathbf{k},p,+\rangle$"))

for i in range(3):
    deltay = 0.25
    deltax = 0.25
    x0,y0 = breite/2+0.5, deltay+dicke
    c.stroke(path.line(x0+i*deltax,y0-deltay, x0+i*deltax+0.2,y0+deltay), [ deco.earrow([deco.stroked([deco.stroked.clear])])])
c.insert(text.text(breite-1.8, 0.5+hoehe/2, r"$|\mathbf{k},p,+\rangle$"))

c.insert(text.text(breite-0.2, -0.3, r"1"))
c.insert(text.text(breite-0.2, hoehe+0.2, r"2"))

c.writePDFfile("roundtrip_pp.pdf")
