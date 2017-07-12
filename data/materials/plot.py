import numpy as np
from pyx import *

hbar = 6.582119514e-16 # eV*s/rad

gradient = color.lineargradient_hsb(color.hsb(1, 1, 0.6), color.hsb(0.4, 1, 0.6))

dalvit = np.loadtxt("GoldDalvit.dat")
mixed  = np.loadtxt("GoldMixed.dat")
#palik  = np.loadtxt("GoldPalik.dat")

xmin, xmax = 1e11, 1e18
x2min, x2max = xmin*hbar, xmax*hbar

g = graph.graphxy(
    width = 8,
    key   = graph.key.key(),
    x     = graph.axis.log(title=r"$\omega$ (rad/s)", min=xmin, max=xmax),
    x2    = graph.axis.log(title=r"$\omega$ (eV)",    min=x2min, max=x2max),
    y     = graph.axis.log(title=r"$\epsilon(i\omega)$", min=1, max=1e8)
)

omegap, gamma = 9/hbar, 0.035/hbar
drude = "y(x) = 1+%g/(x*(x+%g))" % (omegap**2, gamma)

attrs = graph.style.line([gradient, style.linewidth.thick, attr.changelist([style.linestyle.solid, style.linestyle.dashed, style.linestyle.dashdotted])])

g.plot([
    graph.data.points(dalvit, x=1, y=2, title="Palik"),
    graph.data.points(mixed,  x=1, y=2, title="Olmon + Palik"),
    graph.data.function(drude, title="Drude")
], [attrs])

g.writePDFfile()
