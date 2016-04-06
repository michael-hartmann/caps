from __future__ import division
from math import sqrt
from pyx import *

text.set(mode="latex")

attr = [
    deco.earrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
    deco.barrow([style.linewidth.THIn, color.rgb.white, deco.filled([color.rgb.black]), deco.stroked([deco.stroked.clear])]),
]

c = canvas.canvas()

c.stroke(path.circle(0, 0, 1), [style.linewidth.Thin, color.rgb.black, deco.filled([color.cmyk.Periwinkle])])

c.stroke(path.line(-1.4, -2, 1.4, -2), [ style.linewidth.THIck ])

c.writePDFfile("../images/ps.pdf")

c.stroke(path.line(-1.4, -2, 1.8, -2), [ style.linewidth.THIck ])

c.insert(text.text(0.17, 1.35, r"$z$"))
c.stroke(path.line(0, -2.1, 0, 1.5), [ deco.earrow([deco.stroked([deco.stroked.clear])])])


c.insert(text.text(0.1, 0.4, r"$R$"))
c.stroke(path.line(0, 0, 1/sqrt(2), 1/sqrt(2)))

c.stroke(path.line(-1.2, -1, 1.2, -1), [ style.linestyle.dashed ])
c.stroke(path.line(-0.8, 0, 1.8, 0), [ style.linestyle.dashed ])

c.insert(text.text(-0.84, -1.6, r"$L$"))
c.stroke(path.line(-0.54, -2, -0.54, -1), attr)

c.insert(text.text(1.5, -0.9, r"$\mathcal{L}$"))
c.stroke(path.line(1.4, -2, 1.4, 0), attr)

c.writePDFfile("geometry.pdf")
