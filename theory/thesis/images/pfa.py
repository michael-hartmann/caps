from __future__ import division
from math import *
from pyx import *
import numpy as np

phi1 = 255
R = 5
L = 0.8
f = 1.1
text.set(mode="latex")
N = 40

attr = [
    deco.earrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
    deco.barrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
]

c = canvas.canvas()

drawpath = path.path(path.arc(0, R+L, R, phi1, 0))
drawpath.append(path.lineto(R*cos(radians(phi1)), R+L))
drawpath.append(path.lineto(R*cos(radians(phi1)), L+R+R*sin(radians(phi1))))
c.stroke(drawpath, [style.linewidth.Thin, color.rgb.black, deco.filled([color.cmyk.Periwinkle])])

#c.stroke(path.line(R*cos(radians(phi1)), R+L, R, R+L), [ style.linewidth.Thick ])
#c.stroke(path.line(R*cos(radians(phi1)), R+L, R*cos(radians(phi1)), L+R+R*sin(radians(phi1))), [ style.linewidth.Thick ])

c.stroke(path.line(R*cos(radians(phi1)), 0, R*f, 0), [ style.linewidth.THIck ])
c.insert(text.text(4, 5.2, r"$R$"))

c.insert(text.text(5, 4.2, r"$\mathrm{d}x\mathrm{d}y$"))

c.stroke(path.line(0, R+L, R*cos(radians(10)), R+L-R*sin(radians(10))), [ style.linewidth.normal ])

c.stroke(path.line(0, L, R*cos(radians(phi1)), L), [ style.linestyle.dashed ])
c.insert(text.text(-1.3, 0.25, r"$L$"))
c.stroke(path.line(-0.94, L, -0.94, 0), attr)

dx = R/N
i = 0
while True:
    x = i*2*dx
    y = sqrt(R**2-x**2)
    c.stroke(path.line(x-dx, R+L-y, x+dx, R+L-y), [ style.linewidth.Thick ])

    c.stroke(path.line(x-dx, R+L-y, x-dx, 0), [ style.linewidth.Thin, style.linestyle.dashed ])
    #c.stroke(path.line(x-dx, R+L-y, x+dx, R+L-y), [ style.linewidth.Thick ])

    i += 1
    if i*2*dx >= R:
        c.stroke(path.line(x+dx, R+L-y, x+dx, 0), [ style.linewidth.Thin, style.linestyle.dashed ])

        c.stroke(path.line(x, R+L-y, x, 0), attr)
        tbox = text.text(x+0.1, 1.14, r"$d(x,y)$")
        tpath = tbox.bbox().path()
        c.draw(tpath, [deco.filled([color.cmyk.White])])
        c.insert(tbox)
        #c.insert(text.text(x+0.1, 1.4, r"$d(x,y)$"))
        break



c.writePDFfile("pfa.pdf")
