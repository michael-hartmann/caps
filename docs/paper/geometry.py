from math import sin, cos, radians
from pyx import *

R = 1
L = 1/2

text.set(text.LatexRunner)

c = canvas.canvas()
fillstyle = [style.linewidth.thin, color.rgb.black, deco.filled([color.grey(0.8)])]

# sphere
c.stroke(path.circle(0, 0, R), fillstyle)

# plane
c.stroke(path.rect(-2, -R-L, 4, -0.15), fillstyle)

# z-axis
c.insert(text.text(0.17, 1.35, r"$z$"))
c.stroke(path.line(0, -R-L, 0, 1.5), [deco.earrow])

# line and label for radius R (25°)
c.insert(text.text(0.45, -0.05, r"$R$", [text.halign.center]))
phi = radians(-25)
c.stroke(path.line(0, 0, cos(phi), sin(phi)), [deco.earrow, style.linewidth.thick])

# line with arrows and label L
c.insert(text.text(-0.34, -R-L/2, r"$L$", [text.halign.flushcenter, text.valign.middle]))
c.stroke(path.line(0, -R-L, 0, -1), [deco.barrow.small, deco.earrow.small, style.linewidth.thick])

c.writePDFfile()
